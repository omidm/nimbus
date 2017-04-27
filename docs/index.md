# Nimbus Overview

Nimbus is a general purpose cloud computing system specially designed for
computations with short tasks. Nimbus is developed in C++ and the API offers a
data model similar to DraydLINQ and Spark. The key difference between Nimbus
and its counterparts is the novel control plane abstraction called **Execution
Templates**. For jobs with short tasks the runtime overhead becomes comparable
to the task itself. For example for an application with 4ms tasks, Sparks
runtime overhead is about 98%. Execution templates enable Nimbus to schedule
and drive jobs with tasks as short as 100us over large number of workers at
scale with negligible runtime overhead.

Traditional machine learning benchmarks such as logistic regression and k-means
run 40x faster under Nimbus compared to beast available frameworks. In
addition, handling short tasks opens a new class of applications, traditionally
aimed for HPC clusters, into cloud computing world. Nimbus provides a simple
API for physical simulations based on geometry. For example we have ported
PhysBAM, a physics  based simulation library, into Nimbus. Nimbus runs
particle-levelset simulations within 15% of the MPI hand-tuned implementations. 


# Downloading and Building

Nimbus project repository including the core engine source code along with a
handful of applications is available at
[Github](https://github.com/omidm/nimbus).
For complete installation guide please visit [Building Nimbus](building.html).
For a basic installation to run the interactive examples you can issue make at
the Nimbus root:

    $ cd <path-to-nimbus-root>
    $ make


# Running Examples

After you download and build Nimbus you should be able to run simple examples on
your local machine. There are few examples that you can run 


# Where to Go from Here

### Programming Guides:

[Quick Start](quick-start.html): a quick introduction to Nimbus API and how to
write and run your first simple Nimbus application.

[Nimbus Programming Guide](programming-guide.html): detailed overview of Nimbus
API and how to build and link applications against Nimbus core library.


### Deployment Guides:

[Launching Nimbus Cluster:](launching.html) overview of the scripts to launch Nimbus controller
and workers locally or over a cluster of machines. 

[Amazon EC2](ec2.html): scripts that help you launch an elastic cluster of Nimbus nodes
on Amazon EC2. In addition, there are scripts for submitting, monitoring and
managing the application on EC2. 


### Other Documents:

[Monitoring](monitoring.html): how to monitor and profile Nimbus runtime and
application performance.


[Evaluation](evaluation.html): overview of the key benchmark evaluations and
instructions on how to regenerate the results.
