# Multithreading reproducibility test in GitHub Actions CI

Run ISSL-140 example with 1 thread and 4 threads to check that results are bitwise identical.
APCEMM is run for 3h and results are saved every 30 mins to avoid having too many files to compare.

3h seems like enough time for floating point difference to propagate sufficiently to lead to a difference
in outputs without making the CI too long.
