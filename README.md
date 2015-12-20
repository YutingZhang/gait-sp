# Gait-SP: Accelerometer-based gait recognition using signature points

This is the code for the following paper:

Yuting Zhang, Gang Pan, Kui Jia, Minlong Lu, Yueming Wang, Zhaohui Wu, “Accelerometer-based Gait Recognition by Sparse Representation of Signature Points with Clusters”, *IEEE Transactions on Cybernetics*, vol. 45, no. 9, pp. 1864 – 1875, Sept. 2015.


Preparation
===

Compile MinMaxSelection toolbox:

	$ cd bin/MinMaxSelection/MinMaxSelection
	$ matlab
	> minmax_install.m

Compile OMP toolbox:

	$ cd bin/ompbox/private
	$ matlab
	> make

Get extra data, including a MAT version of the [ZJU-GaitAcc dataset](http://www.ytzhang.net/datasets/zju-gaitacc) and some cached experimental results (which can be removed when running from scratch):

	$ ./get_data.sh	


Run
===

	$ cd script
	$ matlab
	> main

After this, you can run the scripts in the current folder to replicate the results in our paper.

Remark
===

This is a beta version. Please feel free to report bugs to me. 
