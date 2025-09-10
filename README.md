# Traffic Assignment Framework with (Customizable) Contraction Hierarchies and Customizable Tree Labeling

This repository contains C++20 code for a traffic assignment framework and optimized 
point-to-point routing algorithms based on state-of-the-art customizable shortest-path algorithms.

The framework and implementation of all shortest-path algorithms except CTL are a fork 
of https://github.com/vbuchhold/routing-framework.
Central ideas for the framework and the usage of CCHs in TA are described in the paper
* [1] Valentin Buchhold, Peter Sanders, and Dorothea Wagner. Real-time Traffic Assignment Using
Engineered Customizable Contraction Hierarchies. ACM Journal of Experimental Algorithmics,
24(2):2.4:1-2.4:28, 2019. [doi:10.1145/3362693](http://dx.doi.org/10.1145/3362693).

The CTL algorithm is based on the paper
* [2] Muhammad Farhan, et al. Customization Meets 2-Hop Labeling: Efficient Routing in Road Networks.
* Proceedings of the VLDB, Volume 18, No. 10, 2025. https://vldb.org/pvldb/volumes/18/paper/Customization%20Meets%202-Hop%20Labeling%3A%20Efficient%20Routing%20in%20Road%20Networks

## License

All files in this repository except the files in the directory `External` are licensed under the MIT
license. External libraries are licensed under their respective licenses. 

## Prerequisites

To build the framework, you need to have some tools and libraries installed. On Debian and its derivatives
(such as Ubuntu) the `apt-get` tool can be used:

```
$ sudo apt-get install build-essential
$ sudo apt-get install cmake
$ sudo apt-get install libboost-all-dev
$ sudo apt-get install sqlite3 libsqlite3-dev
$ sudo apt-get install zlib1g-dev
```

Next, you need to clone the libraries in the `External` subdirectory and build the `RoutingKit` library. To do so,
type the following commands at the top-level directory of the framework:

```
$ git submodule update --init
$ make -C External/RoutingKit lib/libroutingkit.so
```

## Traffic Assignment

### Building Traffic Assignment

Once you installed the packages, type the following commands at the top-level directory of
the framework:

```
$ cmake -S . -B Build/Release -DCMAKE_BUILD_TYPE=Release [-D<option>=<value> ...]
$ cmake --build Build/Release --target AssignTraffic
```

This builds an executable `AssignTraffic` at `Build/Release/Launchers/AssignTraffic`.

There are a number of additional options for routing algorithms that can be passed as CMake parameters:

| Option                          | Possible Values | Default value | Description                                                                                                                                                               |
|---------------------------------|-----------------|---------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `NUM_THREADS`                   | integer         | 1             | Number of threads to use for parallelization of customization and queries in TA. Default is 1.                                                                            |
| `TA_LOG_K`                      | integer         | 0             | Number of simultaneous shortest-path computations using centralized queries (see [1]).                                                                                    |
| `TA_USE_SIMD_SEARCH`            | ON / OFF        | ON            | If set, may use SIMD instructions for centralized shortest-path queries. Only has an effect if `TA_LOG_K` >= 2.                                                           |
| `CTL_THETA`                     | integer         | 0             | Set parameter 'theta', i.e., maximum size of truncated subtrees, of the CTL algorithm (see [2]). Set to 0 for no truncation.                                              | 
| `CTL_SIMD_LOGK`                 | integer         | 0             | Given value L >= 2, the CTL algorithm uses SIMD vector instructions of width 2^L during customization. |
| `CTL_USE_PERFECT_CUSTOMIZATION` | ON / OFF        | OFF           | If set, use perfect customization in CCH on which CTL is based. Default is 'OFF'.                                                                                         |


### Running Traffic Assignment

The `AssignTraffic` executable takes input arguments as outlined at the top of 
the file `Launchers/AssignTraffic.cc`.

The central inputs are a road network in the binary representation of this framework,
and a set of origin-destination (O-D) pairs in a CSV file.

Unfortunately, the input data from the papers mentioned above is not publicly available.
However, we provide scripts to utilize your own input data.
For the road network, consider converting from other formats using the `ConvertGraph` tool in the `RawData` directory.
To map a set of O-D pairs onto a road network in the format of this framework, you can use the `RawData/TransformLocations` tool.

## Point-to-Point Shortest-Path Queries

### Building Point-to-Point Shortest-Path Queries
The framework also provides a launcher for point-to-point shortest-path queries based on the routing algorithms 
implemented. 
Use this launcher to compare the preprocessing, customization, and query performance of the different routing algorithms.
To build the launcher, type the following commands at the top-level directory of the framework:

```
$ cmake -S . -B Build/Release -DCMAKE_BUILD_TYPE=Release [-D<option>=<value> ...]
$ cmake --build Build/Release --target RunP2PAlgo
```

This builds an executable `RunP2PAlgo` at `Build/Release/Launchers/AssignTraffic`.

There are a number of additional options for routing algorithms that can be passed as CMake parameters:

| Option                          | Possible Values | Default value | Description                                                                                                                                                                                                          |
|---------------------------------|-----------------|---------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `NUM_THREADS`                   | integer         | 1             | Number of threads to use for parallelization of customization. (Queries are run sequentially.) Default is 1.                                                                                                         |
| `CTL_THETA`                     | integer         | 0             | Set parameter 'theta', i.e., maximum size of truncated subtrees, of the CTL algorithm (see [2]). Set to 0 for no truncation.                                                                                         | 
| `CTL_SIMD_LOGK`                 | integer         | 0             | Given value L >= 2, the CTL algorithm uses SIMD vector instructions of width 2^L during customization. |
| `CTL_USE_PERFECT_CUSTOMIZATION` | ON / OFF        | OFF           | If set, use perfect customization in CCH on which CTL is based. Default is 'OFF'.                                                                                                                                    |

### Running Preprocessing, Customization and Point-to-Point Shortest-Path Queries

The `RunP2PAlgo` executable takes input arguments as outlined at the top of
the file `Launchers/RunP2PAlgo.cc`.

The central inputs are a road network in the binary representation of this framework,
and, when running queries, a set of origin-destination (O-D) pairs in a CSV file.
The `ConvertGraph` tool in the `RawData` directory can be used to convert a road network from other 
formats to the binary representation of this framework.
For example the road networks from the 9th DIMACS Implementation Challenge on Shortest Paths 
(available at https://www.diag.uniroma1.it/~challenge9/) can easily be converted.
The `GenerateODPairs` tool in the `RawData` directory can be used to generate a set of random O-D pairs in a road network
based on a distribution of O-D distances or Dijkstra ranks.

To evaluate query performance, first use `RunP2PAlgo` to run the preprocessing phase of the respective algorithm and
write the resulting data to a file (see top of file).
Then run `RunP2PAlgo` again with a demand file (`-d` flag) and pass the result of the preprocessing step.

To evaluate customization performance, first run the preprocessing phase as described above.
Then run `RunP2PAlgo` with `CCH-Custom` or `CTL-Custom` for the algorithm parameter.
