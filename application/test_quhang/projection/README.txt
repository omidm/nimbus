Author: Hang Qu <quhang@stanford.edu>
This is a dirty version of PhysBAM projection application.

What this application does is:
(1) Initialize velocity fields.
(2) Calcualte A and b, and other parameters.
(3) Do projection, and solve the pressure.
(4) Use the pressure to update velocity.
(5) Write out velocity for debugging.

This application can run on top of Nimbus, and the result looks good for 2D for two partitions.
It involves using MPI for communication, and uses various job abstraction that cannot be worse.
This application assumes only one region to do projection on, but interfaces to allow for multiple regions are left for future changes.

