# Chase Abram
# Used for testing how to run on Midway

io = open(ENV["SLURM_ARRAY_TASK_ID"]*"test.txt", "w");

write(io, "Test on "*ENV["SLURM_ARRAY_TASK_ID"]);

close(io);

