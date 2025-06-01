# 📦 New Sparse Matrix Storage Format: Contiguous Clustering (CC)

This project introduces **Contiguous Clustering (CC)** — a novel storage format optimized for **diagonally dominant sparse matrices**. It was developed as part of our **Bachelor of Technology thesis** at **XIM University** under the guidance of **Dr. Chandan Misra**.

Traditional formats like CSR, COO, CDS, and JDS each have limitations in space and time efficiency, especially when dealing with matrices that exhibit strong diagonal dominance. Our **CC format** addresses these with smarter clustering and minimal index storage.

---

## 🧠 Motivation

Sparse matrices arise in many real-world applications:
- Scientific simulations
- Social network analysis
- Recommendation systems
- Machine learning and deep learning

However, traditional storage formats often:
- Suffer from excessive memory overhead
- Have inefficient access patterns
- Underperform in SpMV (Sparse Matrix-Vector Multiplication)

---

## 🆕 What is Contiguous Clustering (CC)?

A new storage method that:
- Clusters non-zero values along diagonals
- Stores only **start row/column indices** instead of every (row, col) pair
- Improves cache locality and reduces index redundancy

### Key Components:
- `storeValues[]` – Non-zero matrix values
- `clusterSizes[]` – Number of elements in each diagonal cluster
- `startRowClus[]`, `startColClus[]` – Starting index for each cluster
- `offset[]` – Diagonal offset (col - row)

---

## 📈 Performance Highlights

| Matrix | Format | Space (MB) | Time (s) | GFlops |
|--------|--------|------------|----------|--------|
| mc2dmpi | CDS | 28.0 | 0.12 | 0.03 |
| mc2dmpi | CC  | **12.0** | **0.009** | **0.3** |

- ✅ 30–50% reduction in memory usage
- ✅ Significant improvement in SpMV time
- ✅ Higher GFlops performance across tested datasets

---

## 💻 Implementation

Implemented in **C** with benchmark support for:
- COO Format
- CDS, JDS, PDS Formats
- Our Proposed CC Format

### Compilation
```bash
gcc cc_sparse_matrix.c -o cc_sparse -lm
```
### Run
```bash
./cc_sparse
```


## 📘 Thesis Details
**Title**: New Storage Format for Sparse Matrices
**Institution**: XIM University, Bhubaneswar
**Authors**:

- **Arupa Nanda Swain**

- **Vadali S S Bharadwaja**

- **Satyabhusan Sahu**

- **A Anushruth Reddy**

Guide: Dr. Chandan Misra
Year: 2025

## 🔭 Future Scope
- Extend CC for general sparse matrices (non-diagonal structures)

- Integration with libraries like SciPy, Eigen, PETSc

- GPU-optimized and parallelized versions

- Test on real-world large-scale industrial datasets

## 📜 License
This project is intended for academic and research use only.
© 2025 Arupa Nanda Swain and team. All rights reserved.

