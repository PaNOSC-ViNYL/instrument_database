In this folder, mcstas components should be stored.

Updates from the mcstas repository need to be propagated here from time to time, after checking the impact on the test simulations.


# Testing MCPL input component
Using fish shell syntax.

Create three file for the test and put their names in the list file
```
for s in 0 2 4; mcpltool -e -l2 -s$s /usr/share/mcstas/3.4/data/voutput.mcpl /tmp/voutput-$s.mcpl; echo /tmp/voutput-$s.mcpl.gz >> /tmp/voutput.list; end
mcrun MCPL2hist.instr -d /tmp/MCPL_test -i -n 10
```

