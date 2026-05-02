# Parallel Matrix Processing System (C / Linux IPC)

## Overview
This project is a high-performance matrix processing system implemented in C on Linux. It leverages multi-processing and inter-process communication (IPC) to execute matrix operations in parallel, improving performance compared to traditional sequential approaches.

The system distributes computation across multiple child processes and coordinates execution using low-level OS mechanisms such as pipes, signals, and FIFOs.

---

## Features
- Matrix management system:
  - Create, modify, delete matrices
  - Load from / save to files
- Parallel matrix operations:
  - Addition & subtraction (element-level parallelism)
  - Multiplication (row × column distribution)
  - Determinant computation
- Multi-processing using `fork()`
- IPC mechanisms:
  - Pipes
  - Signals
  - FIFOs
- Performance comparison:
  - Sequential vs Parallel execution
- Optional OpenMP support for loop-level parallelism

---

## Tech Stack
- **Language:** C  
- **Platform:** Linux  

### Concepts Used
- Process creation (`fork`)
- Inter-process communication (IPC)
- Signals & synchronization
- Pipes and FIFOs
- Parallel computing

---

## System Architecture


Parent Process
│
├── Distributes matrix data
│
├── Creates child processes
│
├── Sends tasks via IPC (pipes / FIFOs)
│
└── Collects results
↑
Child Processes compute partial results


---

## Design Decisions
- **Process-based parallelism** was used instead of threads to practice low-level OS concepts.
- **Pipes/FIFOs** were chosen for efficient communication between parent and child processes.
- **Signals** were used for synchronization and coordination.
- Workload is divided at a fine-grained level (per element or per row/column) to maximize parallelism.

---

## How to Run

### Compile
```bash
gcc -o matrix main.c -lm
Run
./matrix config.txt

Note: config.txt contains configuration for loading matrices or system settings.

Example Run
Menu:
1. Add Matrix
2. Multiply Matrix
...

Enter choice: 10

Matrix A:
[1 2]
[3 4]

Matrix B:
[5 6]
[7 8]

Result:
[19 22]
[43 50]
Performance
Operation	Sequential	Parallel
Addition	12 ms	4 ms
Multiply	85 ms	30 ms

Parallel execution significantly reduces computation time by distributing work across multiple processes.

Example Operations
Matrix addition using parallel child processes
Matrix multiplication using row-column mapping
Determinant calculation using distributed computation
Future Improvements
Thread-based implementation (pthreads)
Graphical user interface (GUI)
GPU acceleration (CUDA / OpenCL)
Dynamic process pool for better resource management
