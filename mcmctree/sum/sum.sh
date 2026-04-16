export PATH=/data/software/paml4.9j/src:$PATH

head -n1 ../gtr_approx01/mcmc.txt > mcmc.txt
tail -n800000 ../gtr_approx01/mcmc.txt >> mcmc.txt 
tail -n800000 ../gtr_approx02/mcmc.txt >> mcmc.txt
tail -n800000 ../gtr_approx03/mcmc.txt >> mcmc.txt
tail -n800000 ../gtr_approx04/mcmc.txt >> mcmc.txt

mcmctree mcmctree.ctl > screenlog.log 2>&1
