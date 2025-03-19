# 🛠️ **Generalized Dual Simplex Solver for LPP**

This project implements a **generalized solver** for Linear Programming Problems (LPP) using the **Dual Simplex Method** in Python. It standardizes the given problem, handles **unrestricted variables**, and displays the **primal, dual, and standard dual forms** before solving.

---

### 🔥 **Features**
- ✅ **Automatic Standardization:** Converts the given LPP to its standard form.  
- 🔁 **Dual Conversion:** Automatically generates the dual problem.  
- 🔓 **Unrestricted Variables:** Handles variables with no sign restrictions using the `x' - x''` transformation.  
- 📊 **Dual Simplex Method:** Solves the problem iteratively and displays the tableau.  
- 🎯 **Optimal Solution Display:** Prints the optimal solution or declares infeasibility.  

---

### 🚀 **How to Run**
1. **Clone the repository:**  
```bash
git clone https://github.com/yourusername/generalized-dual-simplex-lpp.git
cd generalized-dual-simplex-lpp
```
2. **Install dependencies:**  
```bash
pip install numpy pandas
```
3. **Run the solver:**  
```bash
python main.py
```

---

### 📄 **Usage Example**
```  
Is the input in primal form? (yes/no): yes  
Enter max or min: max  
Enter number of decision variables: 3  
Enter number of constraints: 2  
Enter objective coefficients: 3 2 4  
Constraint 1 coefficients: 2 3 1  
Constraint 2 coefficients: 4 1 2  
Enter signs: <= <=  
Enter RHS values: 5 8  
Enter unrestricted variables (e.g., 1 2): 1  
```
**Output:**  
```
Optimal Solution:  
X1 = 1.25  
X2 = 0  
X3 = 0  
Optimal Z = 3.75  
```

---

### 🛠️ **File Structure**
```
/generalized-dual-simplex-lpp  
 ├── main.py                # Main execution script  
 ├── dual_simplex.py        # Dual simplex algorithm  
 ├── standardize.py         # Primal and dual standardization  
 ├── utils.py               # Helper functions (printing, tableau)  
 ├── README.md              # Project documentation  
 ├── requirements.txt       # Dependencies  
 ├── LICENSE                # License file  
 └── .gitignore             # Git ignore file  
```

---

### 📚 **Contributing**
Feel free to fork the repository and submit pull requests with improvements or bug fixes.

---

### ⚙️ **License**
[MIT License](LICENSE)
