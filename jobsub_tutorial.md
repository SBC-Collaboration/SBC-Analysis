# Jobsub Setup

This is required to setup batch jobs on the virtual machine.

#### Kerberos Ticket Forwarding
If you are using PuTTy or different kerberos ticket manager that isn't bash, skip to the end of this section.

The batch job manager runs on an internal server at fermilab that requires kerberos authentication. From your host computer, allow your kerberos ticket to be forwarded to the virtual machine when it is initialized. This is done by adding the -f option when instantiating your ticket and the -K option when sshing into the virtual machine. You can then check your active tickets by running the klist command.
```
$ kinit -f your_username@FNAL.GOV
$ ssh -K your_username@couppgpvm01.fnal.gov
$ klist
```

If you are using PuTTy or another external kerberos ticket manager, there may be an option to forward your ticket, otherwise you will just have to create a second kerberos ticket once you've logged onto the virtual machine. I

#### Jobsob Access
[jobsub](https://cdcvs.fnal.gov/redmine/projects/jobsub/wiki/Using_the_Client) is an external package that allows users to utilize much more RAM/CPU time/Disk space than they normally would for completing scientific tasks. To gain access to the jobsub suite, you will have to run this every time you open a shell or simply add it to a startup script.
```
$ . /grid/fermiapp/products/common/etc/setups.sh
$ setup jobsub_client
```

#### Create and Submit a Test Job
```
$ nano my_first_jobsub.sh
echo This is my first jobsub script
sleep 1m # <-- Gives us enough time to check status

$ jobsub_submit --group=coupp file://my_first_jobsub.sh --expected-lifetime=short
```
After submitting your job, you should receive a confirmation message along with a jobid to fetch your job/check the status. 
#### Checking Job Status
Let's assume your jobid was `1032512.0@jobsub01.fnal.gov`.

* Check the status of your job with `jobsub_q`
```
$ jobsub_q --group=coupp --jobid=1032512.0@jobsub01.fnal.gov
``` 
* Collect your jobs output with `jobsub_fetchlog`. For large batch jobs, there may be thousands of output files, so it's wise to create an output directory.
```
$ mkdir joboutput
$ cd joboutput
$ jobsub_fetchlog --group=coupp --jobid=1032512.0@jobsub01.fnal.gov
$ tar -xzvf TheFileNameFrom^ .
```
