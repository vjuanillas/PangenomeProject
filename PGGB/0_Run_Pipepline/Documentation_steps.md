# **Documentation on running the PGGB Pangenome on the HPZ server**


# **Steps in running:**

1. ## **Prepare the environment:** 

1. Prepare the working environment. Python version is pinned on 3.10 because PGGB 0.7.2 requires python-igraph 0.11.5, which only works with Python 3.8–3.12.  
   

| $ sh 000\_prepare\_environment.sh |
| :---- |

	

2. Activate the pggb\_pangenome virtual environment

| $ conda activate pggb\_pangenome |
| :---- |

3. Install required programs using the 001\_prepare\_environment.sh script

	

| $ sh 001\_prepare\_environment.sh |
| :---- |

2. ## **Prepare the assembly sequences**

| $ sh 00\_prepare\_sequences.sh |
| :---- |

3. ## **Run PGGB**

   

| $ sh 01\_run\_pggb.sh |
| :---- |

4. ## **Characterize the pangenome**

| $ sh 02\_charactetize.sh |
| :---- |

5. 
