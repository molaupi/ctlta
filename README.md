# Routing Framework

## License

All files in this repository except the files in the directory `External` are licensed under the MIT
license. External libraries are licensed under their respective licenses. Note that the compiled
programs `CalculateDemand` and `ComputeUnionBoundary` use libraries that are released under the GNU
GPLv3, and thus the compiled programs `CalculateDemand` and `ComputeUnionBoundary` have to be under
the GNU GPLv3.

## Prerequisites

To build the framework, you need to have some tools and libraries installed. On Debian and its derivatives
(such as Ubuntu) the `apt-get` tool can be used:

```
$ sudo apt-get install build-essential
$ sudo apt-get install cmake
$ sudo apt-get install python3 python3-pip; pip3 install -r python_requirements.txt
$ sudo apt-get install sqlite3 libsqlite3-dev
$ sudo apt-get install zlib1g-dev
$ sudo apt-get install libproj-dev
```

Next, you need to clone the libraries in the `External` subdirectory and build the `RoutingKit` library. To do so,
type the following commands at the top-level directory of the framework:

```
$ git submodule update --init
$ make -C External/RoutingKit lib/libroutingkit.so
```

## Building

Once you installed the packages, type the following commands at the top-level directory of
the framework:

```
$ cmake -S . -B Build/Release_<theta> -DCMAKE_BUILD_TYPE=Release -DTTL_THETA=<theta>
$ cmake --build Build/Release_<theta> --target AssignTraffic
```

This builds an executable `AssignTraffic` using a theta value of
`<theta>` at `Build/Release_<theta>/Launchers/AssignTraffic`.
Replace `<theta>` with the desired theta value. 
Run the two commands multiple times for multiple `AssignTraffic` executables with 
different theta values in different directories. 

## Running Traffic Assignment

The `AssignTraffic` executable takes input arguments as outlined at the top of 
the file `Launchers/AssignTraffic.cc`.
For example, for `Stuttgart-morn-undirected`, use:

```
$ Build/Release_<theta>/Launchers/AssignTraffic -v -p 1 -n 50 -a TTLSA \
    -g <path>/Inputs/Graphs/Visum_Stuttgart_undirected.gr.bin
    -d <path>/Inputs/ODPairs/Visum_Stuttgart_morn_undirected.csv
    -flow <path>/Outputs/Stuttgart_morn_undirected/flow_TTLSA_TH<theta> \
    -dist <path>/Outputs/Stuttgart_morn_undirected/dist_TTLSA_TH<theta> \
    -stat <path>/Outputs/Stuttgart_morn_undirected/stat_TTLSA_TH<theta>
```

The argument `-v` makes the executable print out statistics for every 
iteration to standard out.

The argument `-p` sets the observation period in hours. It should be set to 1 
for `morn` and `even`, to 24 for `day` and to 168 for `week`.

The argument `-n` specifies the number of iterations to run.

The argument `-a` specifies the algorithm. TTLSA stands for TTL standalone, 
which is your implementation of CTL. It is run with the theta value specified 
earlier when building. 
If you want to run TA with your CCH implementation use `-a TTLSACCH`. 
The theta value doesn't matter for this, obviously. 

The argument `-g` is a path to the road network in the binary representation 
of this framework. You should already have received the Stuttgart graphs in this 
format. If you need to convert from DIMACS to this format, take a look at 
`RawData/ConvertGraph.cc` (or ask me).

The `-d` argument specifies a path to the O-D pairs as a CSV file with a column called 
origin and a column called destination containing the vertex IDs of the 
origins and destinations, respectively. You should have received the Stuttgart 
demand sets.

The `-flow`, `-dist`, and `-stat` arguments specify paths to the three output 
files in CSV format. 
The flow output file contains the flow on each edge after the first 
and after the last iteration.
The dist output file contains the travel cost for each O-D pair in the first
and in the last iteration.
The stat output file contains statistics on the runtime of each iteration 
in milliseconds. You'll be interested in the `query_time` and 
`customization_time` columns.