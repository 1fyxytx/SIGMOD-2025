# SIGMOD-2025
The source code of ParaEnum is QualityPath.cpp. To run this program, you need to do the following steps:
1. Download the datasets from "http://konect.cc/networks/", "https://networkrepository.com/index.php", and "https://law.di.unimi.it/", which have the same format as test.graph.
2. Run the done.sh file.

The explanations of some parameters are listed as follows:
1. weig: the maximal weight value in the dataset.
2. threads: the number of computing cores.
3. TaskCnt: the number of query tasks.
4. MaxLoad: the pre-defined maximal skyline paths of the query tasks that are executed in a single core.
