🔹 Project Title

Parallel Matrix Processing System (C / Linux IPC)

🔹 Description

This project is a high-performance matrix processing system implemented in C using Linux multi-processing and inter-process communication (IPC) mechanisms.

The system performs various matrix operations (addition, subtraction, multiplication, determinant, etc.) by distributing computations across multiple child processes to achieve parallel execution.

🔹 Features

Matrix management system (create, modify, delete, load/save)
Parallel matrix addition and subtraction using child processes
Parallel matrix multiplication (row × column distribution)
Determinant computation using process-based parallelism
Inter-process communication using:
Pipes
Signals
FIFOs
Performance comparison between:
Sequential execution
Parallel execution
Optional OpenMP optimization for loop parallelism

🔹 Tech Stack

Language: C
Platform: Linux
Concepts:
fork()
signals
pipes
FIFOs
process synchronization
parallel computing

🔹 System Design

Each matrix operation is divided into smaller tasks
Tasks are distributed across multiple child processes
Parent process coordinates:
data distribution
result collection
IPC mechanisms are used for communication between processes

🔹 How to Run

gcc -o matrix main.c -lm
./matrix config.txt



🔹 Example Operations

Add two matrices using parallel processes
Multiply matrices using row-column mapping
Compute determinant using distributed 
Example:
Matrix A:
[1 2]
[3 4]

Matrix B:
[5 6]
[7 8]

Result:
[19 22]
[43 50]


## Performance

| Operation | Sequential | Parallel |
|----------|-----------|----------|
| Addition | 12 ms     | 4 ms     |
| Multiply | 85 ms     | 30 ms    |


Parent Process
   ↓
Distributes tasks
   ↓
Child Processes compute elements
   ↓
Results sent via pipes 

🔹 Future Improvements

Thread-based implementation (pthreads)

GUI interface

GPU acceleration (CUDA/OpenCL)

Dynamic process pool optimization
