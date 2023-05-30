# 3D Contrast Enhanced Ultrasound (CEUS) Analysis GUI

## Overview

This program is a CEUS analyis tool which allows user to input a string of 3D B-Mode images, draw a volume of interest, and view the resulting time intensity curve (TIC).

Next, the user can choose the point to start analysis on the TIC (at point t0), modify the TIC to remove noise from the data, and fit a lognormal curve to the resulting TIC.

### Initial TIC Image

![Initial TIC Image](images/initTICImage.png "Initial TIC Image")

### Modifying TIC

![Editing TIC Image](images/midTICImage.png "Modifying TIC Image")

### Final TIC Image

![Final TIC Image](images/finalTICImage.png "Final TIC Image")

Finally, using this lognormal curve fitting, the program computes the normalized area under the curve (AUC), normalized peak value (PV), normalized time to peak (TP), normalized mean transit time (MTT), and normalizing value (TMPPV).

### Main GUI

![Main GUI Image](images/imageGUI.png "Main GUI Image")

## Dependencies

* [Python Version 3.9.13](https://www.python.org/downloads/release/python-3913/)
* [Git](https://git-scm.com/downloads)

## Building

### Mac/Linux

```shell
git clone https://github.com/davidspector67/3dCEUSAnalysis.git
cd 3dCEUSAnalysis
chmod +x init.sh
chmod +x run.sh
./init.sh
```

### Windows

```shell
git clone https://github.com/davidspector67/3dCEUSAnalysis.git
cd 3dCEUSAnalysis
pip install virtualenv
virtualenv venv
.\venv\Scripts\activate.bat
pip install -r .\pyPackages.txt
```

## Running

### Mac/Linux

```shell
./run.sh
```

### Windows

```shell
.\run.sh
```

```shell
./hash-table-tester -t [NUM_CORES] -s [NUM_ENTRIES]
# Generation: 36,996 usec
# Hash table base: 84,249 usec
#   - 0 missing
# Hash table v1: 191,023 usec
#   - 0 missing
# Hash table v2: 22,691 usec
#   - 0 missing
```

Here, this lab can be run using the above command format where `NUM_CORES` is a positive integer denoting the number
of threads to be used and `NUM_ENTRIES` is a positive integer denoting the number of hash table entries to be used by each thread to test each "add entry" implementation.

## First Implementation

In the `hash_table_v1_add_entry` function, I added the `mutex` at line 95. I blocked the entire `hash_table_v1_add_entry` function as a critical section, guaranteeing correctness because all operations are protected which prevents any potential race condition between competing threads.

### Performance

```shell
./hash-table-tester -t 4 -s 50000
# Generation: 28,666 usec
# Hash table base: 78,995 usec
#   - 0 missing
# Hash table v1: 200,597 usec
#   - 0 missing
```

Version 1 is slower than the base version. Since the critical section is the entirety of the `hash_table_v1_add_entry` function, only one function call can be handled at a time by the CPU. Thus, Version 1 must complete the same amount of instructions sequentially as the base version, except Version 1 must also create multiple threads and initialize a mutex. This causes additional overhead in the program and thus, Version 1 takes more time to complete than the base version. Specifically, thread creation means that a new stack needs to be initialized for each thread and some register values must also be stored and loaded during each thread switch. Since we never use the concurrent properties of threads, this process can be thought of as added overhead to our program without any added benefit to go with it.

I initalized the mutex I used globally on line 12 and destroyed it in the `hash_table_v1_destroy` function.

## Second Implementation

In the `hash_table_v2_add_entry` function, I used one mutex per hash table entry. Each mutex protects the linked list corresponding to each hash table entry, and ensures that new data entries to the linked list are added atomically.

First, I added these mutexes into the `hash_table_entry` data structure and initialized them in the `hash_table_v2_create` function.

Next, within the `hash_table_v2_add_entry` function, my single critical section contained the attempted access of the inputted list entry, the update of an existing list entry's value, and the addition of a new data entry into the linked list in the current hash table entry. The attemped access of the inputted list entry via the `get_list_entry` function is in the critical section, because it is a dependent read for our write to the shared hash table. Similarly, the if statement on line 103 is also a dependent read for our write to the shared hash table. Next, both the data update at line 104 and the linked list insert on line 116 are shared writes which rely on the same dependent reads, so they're also included in the critical section as well.

Next, although the list_entry allocation and key and value assignments (lines 112-114) are thread-safe since they deal with a local variables, they are still included in the critical section. This is because as previously mentioned, the linked list insertion call at line 116 must be in our one critical section, and I made the choice to keep this code in the critical section which blocks infrequently rather than to always allocate memory before the critical section for a data entry that may not be used.

Finally, all used mutexes were destroyed in the `hash_table_v2_destroy` function.

No more locks are required, because all possible data races occur at the same hash table entry.

### Performance

```shell
./hash-table-tester -t 4 -s 50000
# Generation: 36,996 usec
# Hash table base: 84,249 usec
#   - 0 missing
# Hash table v1: 191,023 usec
#   - 0 missing
# Hash table v2: 22,691 usec
#   - 0 missing
```

This time the speed up is approximately 3.7 times faster, similar than the physical core number 4. This is because all parts of the `hash_table_v2_add_entry` function outside of the previously described critical section are run entirely in parallel, saving us time proportional to the number of cores we have. Also, when multiple threads aren't trying to add a data entry to the same hash table entry simultaneously (which occurs the majority of the time), they can all run their critical sections in parallel, because none of them share any mutexes.

Thus, the entirety of the `hash_table_v2_add_entry` function runs in parallel with the exception of the case where two threads are indeed trying to add a data entry to the same hash table entry simulataneously. In this case, only one thread can run its critical section at once since only one mutex is shared between the competing threads. However, since this case is rare, we have that most of the time, all of our cores can be running the entirety of `hash_table_v2_add_entry` simultaneously. Since we have 4 cores, this comes out to roughly a 4 times sped up. Our actual speed up of 3.7 times is due in part to the occasional thread(s) that must wait to use its mutex in the previously described rare case. The overhead of managing multiple threads and mutexes also contributes to our actual speed up being less than 4.

## Cleaning up

```shell
make clean
```
