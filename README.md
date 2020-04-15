
# Nimbus Project

Nimbus is a framework for fast cloud computing applications. It is developed
by a group of researchers at Stanford University. Nimbus supports jobs with
short tasks and high task rate. The range of applications includes machine
learning, graph processing and graphical simulations. Nimbus's control plane
implements a novel abstraction, called **execution termplates** that allows
Nimbus to scale out and provide orders of magnitude higher task throughput
compared to counterpart frameworks. Nimbus introduces a simple geometric based
interface and a translation layer for automatic distribution of graphical
simulations in the cloud. For more information please refer to [Nimbus
website](http://nimbus.stanford.edu).

Followings are concise instructions on how to build and use Nimbus. For the
complete instructions please refer to the `docs` folder or access an online
version of the [Nimbus documentation](https://omidm.github.io/nimbus/). 


## Building Nimbus

Nimbus is built using `make` and requires a `gcc` compiler. We have tested
Nimbus with `gcc` versions 4.5, 4.8 and 5.4, on Ubuntu 12.04, 14.04 and 16.04.
Also, you might be able to build Nimbus on OSX (although it is not officially
supported). To build Nimbus, first install `make` and `gcc` on your system (you
most probably already have them). For example, on Ubuntu:

    $ sudo apt-get install gcc
    $ sudo apt-get install g++
    $ sudo apt-get install make

Then issue make at Nimbus root directory: 

    $ make

It builds the nimbus core library and external library dependencies along with
a group of applications. Nimbus library uses `boost`, `protobuf`, and `leveldb` header
files and libraries. However, you do not need to install them
separately. They are locally installed as part of Nimbus built, automatically.

The build might take few minutes. You could use `-j <number of threads>` option
with `make` for faster turn around. After the build is successful, you can test
basic functionalities by running simple end-to-end applications. For example,
at the root, issue:
  
    $ make test-basics

For complete building instructions and how to build and run graphical
simulations with Nimbus please refer to the documentations,
at [Building Nimbus](https://omidm.github.io/nimbus/building.html).



## Running Examples

You can run basic, one controller one worker, application examples with a
single script.  From the root directory call `./scripts/run-examples.sh`. For
example, you can run logistic regression application as follows:

    $ ./scripts/run-examples.sh lr 

To see other available examples through this script use the `-h` options:
  
    $ ./scripts/run-examples.sh -h


Also, You can run more complex examples easily with one controller against
multiple workers. Use the scripts to launch controller and workers from the
root directory. To launch the controller against two worker, and split the
application domain along x-axis, issue:

    $ ./scripts/start-controller.sh -w 2 --split 2 1 1

This script starts the controller locally and redirects the controller
`stdout/stderr` to files in `logs/` folder. If you want to run the controller in
the foreground and flush the logs in to current console use the `--fg` option.
There are numerous options that you can use for the controller. For the entire
list, issue:

    $ ./scripts/start-controller.sh -h

Now you need to launch the workers to run against the controller. To launch the
two workers to run for example the k-means application, issue:

    $ ./scripts/start-workers.sh 2 -l applications/ml/k_means/libk_means.so

This script launches two workers, and redirects their `stdout/stderr` to files
in  `logs/` folder. To print the logs in the console use the `--fg` option; in
case of multiple workers it opens a terminal for each worker. Also, workers
have multiple options, to see the complete list:

    $ ./scripts/start-workers.sh -h
  

