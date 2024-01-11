for f in *; do echo mv $f `echo $f | sed 's/_R_.*/.amplicon.cov.xls/'`; done
### line line is for echoing the result, remove first echo after do to make it work
