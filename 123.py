import numpy as np
import pandas as pd

def print_primaldual(c, A, s, b, opt_type, unrestricted_vars=None):
    num_constraints, num_vars = len(b), len(c)

    # Use original variable names (e.g., X1, X2, X3) for input display
    X = [f"X{i+1}" for i in range(num_vars)]
    
    print("Maximize") if opt_type == 'max' else print("Minimize")
    
    print("\t", end=" ")
    for i in range(num_vars-1):
        print(c[i], X[i], end=" ")
        print("+", end=" ")
    print(c[num_vars-1], X[num_vars-1]) 
    
    print("\nSubject to:")
    
    for i in range(num_constraints):
        print("\t", end=" ")
        for j in range(num_vars-1):
            print(A[i][j], X[j], end=" ")
            print("+", end=" ")
        print(A[i][num_vars-1], X[num_vars-1], end=" ")
        print(s[i], end=" ")
        print(b[i])
        
    print("\t", end=" ")    
    for i in range(num_vars):
        if unrestricted_vars and (i+1) in unrestricted_vars:
            print(X[i], "is unrestricted", end=", " if i < num_vars-1 else "")
        else:
            print(X[i], ">= 0", end=", " if i < num_vars-1 else "")
    print()         
                
    print("\n" + "="*100 + "\n")


def replace_resunres(c, A, unrestricted_vars, num_variables):
    # Replace unrestricted variables with x' - x''
    new_c = []
    new_A = []
    for i in range(num_variables):
        if i+1 in unrestricted_vars:
            # Add x_i' and x_i'' to the objective function
            new_c.append(c[i])
            new_c.append(-c[i])
        else:
            # Add x_i as is
            new_c.append(c[i])
    
    # Convert new_c to a numpy array
    new_c = np.array(new_c)
    
    # Adjust the constraint matrix A
    for row in A:
        new_row = []
        for j in range(num_variables):
            if j+1 in unrestricted_vars:
                # Add coefficients for x_j' and x_j''
                new_row.append(row[j])
                new_row.append(-row[j])
            else:
                # Add coefficient for x_j as is
                new_row.append(row[j])
        new_A.append(new_row)
    
    # Convert new_A to a numpy array
    new_A = np.array(new_A)
        
    return new_c, new_A

    
def standardize_primal(c, A, s, b, opt_type):
    num_constraints, num_vars = len(b), len(c)

    count = 0    
    for i in range(len(s)): 
        if s[i] == "=":
             count = count+1

    Ap = np.zeros((num_constraints + count, num_vars))
    bp = np.zeros(num_constraints + count)
    sp = [""]*(num_constraints + count)
    
    j = 0
    for i in range(len(s)):
        Ap[j, :] = A[i, :].copy()
        bp[j] = b[i]
        sp[j] = "<="
        if s[i] == "=":
            j = j+1
            Ap[j, :] = -Ap[j-1, :]
            bp[j] = -bp[j-1]
            sp[j] = "<=" 
        elif s[i] == ">=":
            Ap[j] = -Ap[j]
            bp[j] = -bp[j]    
        j = j+1
    
    """Ensures the primal LPP is in max form with ≤ constraints."""
    if opt_type == "min":
        c = -c  # Convert min to max

    return c, Ap, sp, bp, "max"

def convert_to_dual(cp, Ap, sp, bp, opt_type):
    """Converts the standard primal LPP to its dual and ensures correct standardization."""
    num_variables = len(cp)
    
    dual_s = [">="]*(num_variables)
    
    dual_c = np.copy(bp)
    dual_A = np.transpose(Ap)
    dual_b = np.copy(cp)
    
    dual_opt_type = "min"
    return dual_c, dual_A, dual_s, dual_b, dual_opt_type

def standardize_dual(c, A, s, b, opt_type):
    """Ensures the dual LPP is in max form with ≤ constraints."""
    if opt_type == "min":
        c = -c
        A = -A
        b = -b
        for i in range(len(s)):
            s[i] = "<="
        
    return c, A, s, b, "max"


def print_tableau(tableau, basic_vars, C, cost, ratios, pivot_row, pivot_col):
    """Prints the tableau in a structured format."""
    num_constraints, num_vars = tableau.shape[0] - 1, len(C)
    
    # Column names
    columns = ["CB"] + [f"X{i+1}" for i in range(num_vars)] + [f"S{i+1}" for i in range(num_constraints)] + ["RHS"]
    
    df = pd.DataFrame(tableau, columns=columns)
    df.loc[-1] = cost
    df.index = df.index + 1
    df = df.sort_index()

    # Theta (only for non-initial tableaus)
    
    theta = [""] * (num_vars + num_constraints + 2)  
    if len(ratios) != 0:
        for i in range(len(ratios)-2):
            if ratios[i+1] != np.inf:
                theta[i+1] = ratios[i+1] 
                    
    df.loc[num_constraints+2] = theta

    # Pivot column marker
    pivot_col_symb = [""] * (num_vars + num_constraints + 2)
    if pivot_col != -1:
        pivot_col_symb[pivot_col] = '|'    
    df.loc[num_constraints+3] = pivot_col_symb

    # Basic variable column
    basic_vars_col = [""] * (num_constraints + 4)
    basic_vars_col[0] = 'C_j'
    for i, bv in enumerate(basic_vars):
        basic_vars_col[i+1] = bv
    basic_vars_col[num_constraints+1] = 'C_j-Z_j'
    if len(ratios) != 0:
        basic_vars_col[num_constraints+2] = 'theta'    

    df.insert(0, "Basic", basic_vars_col)

    # Pivot row marker
    pivot_row_symb = [""] * (num_constraints + 4)
    if pivot_row != -1:
        pivot_row_symb[pivot_row+1] = '-->'
    df.insert(num_vars + num_constraints + 3, " ", pivot_row_symb)

    print(df.to_string(index=False))
    print("\n" + "="*100 + "\n")


def dual_simplex(c, A, b, opt_type):
    """Runs the dual simplex method."""
    num_variables = len(c)
    num_constraints = len(b)
    
    # Convert min to max if necessary
    if opt_type == "min":
        c = -c

    cost = [0.0] * (num_variables + num_constraints + 2)
    cost[0] = ""
    for i in range(num_variables):
        cost[i+1] = c[i]
    cost[num_variables+num_constraints+1] = ""        

    tableau = np.zeros((num_constraints + 1, num_variables + num_constraints + 2))
    tableau[:-1, 1:num_variables+1] = A
    tableau[:-1, num_variables+1:num_variables+num_constraints+1] = np.eye(num_constraints)
    tableau[:-1, -1] = b
    tableau[-1, 1:num_variables+1] = c

    basic_vars = [f"S{i+1}" for i in range(num_constraints)]
    print("\nInitial Tableau:")
    print_tableau(tableau, basic_vars, c, cost, [], -1, -1)

    while np.any(tableau[:-1, -1] < 0):  # Check if any RHS is negative
        pivot_row = np.argmin(tableau[:-1, -1])
        if tableau[pivot_row, -1] >= 0:
            break  

        ratios = []
        for j in range(tableau.shape[1] - 1):
            if tableau[pivot_row, j] < 0:
                ratios.append(abs(tableau[-1, j] / tableau[pivot_row, j]))
            else:
                ratios.append(np.inf)
        
        ratios[0] = np.inf
        pivot_col = np.argmin(ratios)

        if np.all(np.array(ratios) == np.inf):
            print("No feasible solution exists!")
            return "NFSE"

        print("\nPivot Selection:")
        print_tableau(tableau, basic_vars, c, cost, ratios, pivot_row, pivot_col)

        pivot_element = tableau[pivot_row, pivot_col]
        tableau[pivot_row] /= pivot_element

        for i in range(tableau.shape[0]):
            if i != pivot_row:
                tableau[i] -= tableau[i, pivot_col] * tableau[pivot_row]

        tableau[pivot_row, 0] = 0

        if pivot_col <= num_variables:
            tableau[pivot_row, 0] = cost[pivot_col]
           
        if pivot_col <= num_variables:
            basic_vars[pivot_row] = f"X{pivot_col}"
        else:
            basic_vars[pivot_row] = f"S{pivot_col-num_variables}"   

        print("\nUpdated Tableau:")
        print_tableau(tableau, basic_vars, c, cost, [], -1, -1)

    solution = np.zeros(num_variables)
    for i in range(num_variables):
        col = tableau[:, i+1]
        if np.count_nonzero(col[:-1]) == 1 and np.sum(col[:-1]) == 1:
            row_index = np.where(col[:-1] == 1)[0][0]
            solution[i] = tableau[row_index, -1]

    optimal_value = -tableau[-1, -1]
    if opt_type == "min":
        optimal_value = -optimal_value  

    print("\nOptimal Dual Solution:")
    for i in range(num_variables):
        print(f"x{i+1} = {solution[i]:.3f}")
    print(f"Optimal Dual Z = {optimal_value:.3f}")

    return optimal_value

def main():
    form = input("Is the input in primal form? (yes/no): ").strip().lower()
    opt_type = input("Enter max or min: ").strip().lower()
    num_variables = int(input("Enter number of decision variables: "))
    num_constraints = int(input("Enter number of constraints: "))
    
    c = np.array(list(map(float, input("Enter objective coefficients: ").split())), dtype=float)
    A = np.array([list(map(float, input(f"Constraint {i+1} coefficients: ").split())) for i in range(num_constraints)], dtype=float)
    s = np.array(list(map(str, input("Enter signs: ").split())), dtype=str)
    b = np.array(list(map(float, input("Enter RHS values: ").split())), dtype=float)
    
    unrestricted_vars = list(map(int, input("Enter unrestricted variables (e.g., 1 2): ").split()))
    
    print("\n" + "="*100 + "\n")
    print("User Input:\n")
    print_primaldual(c, A, s, b, opt_type, unrestricted_vars)

    if form == "yes":
        # Replace unrestricted variables with x' - x''
        new_c, new_A = replace_resunres(c, A, unrestricted_vars, num_variables)
        
        # Standardize the primal problem
        cp, Ap, sp, bp, opt_type1 = standardize_primal(new_c, new_A, s, b, opt_type)
        
        print("\n" + "="*100 + "\n")
        print("Standard Primal Form:\n")
        print_primaldual(cp, Ap, sp, bp, opt_type1)
        
        # Convert to dual
        cd, Ad, sd, bd, opt_type1 = convert_to_dual(cp, Ap, sp, bp, opt_type1)
         
        print("\n" + "="*100 + "\n") 
        print("Dual Form:\n") 
        print_primaldual(cd, Ad, sd, bd, opt_type1)
 
        # Standardize the dual problem
        cd, Ad, sd, bd, opt_type1 = standardize_dual(cd, Ad, sd, bd, opt_type1)

        print("\n" + "="*100 + "\n") 
        print("Standard Dual Form:\n") 
        print_primaldual(cd, Ad, sd, bd, opt_type1)
        
        # Apply dual simplex if necessary
        
        if np.any(bd < 0):
            print("\nApplying Dual Simplex:")
            optimal_value = dual_simplex(cd, Ad, bd, opt_type1)
            
            if optimal_value != "NFSE":
                if opt_type == "max":
                    optimal_value = -optimal_value
                print(f"Optimal Primal Z = {optimal_value:.3f} or Infeasible")            
        else:
            print("\nAll bi values are positive. No need to apply Dual Simplex. No solution.")
    
    else:
        # Input is already in dual form
        # Standardize the dual problem
        cd, Ad, sd, bd, opt_type1 = standardize_dual(c, A, s, b, opt_type)

        print("\n" + "="*100 + "\n") 
        print("Standard Dual Form:\n") 
        print_primaldual(cd, Ad, sd, bd, opt_type1)
        
        # Apply dual simplex if necessary
        if np.any(bd < 0):
            print("\nApplying Dual Simplex:")
            optimal_value = dual_simplex(cd, Ad, bd, opt_type1)
            
            if optimal_value != "NFSE":
                if opt_type == "min":
                    optimal_value = -optimal_value
                print(f"Optimal Standard Dual Z = {optimal_value:.3f}")            
        else:
            print("\nAll bi values are positive. No need to apply Dual Simplex. No solution.")
    
if __name__ == "__main__":
    main()