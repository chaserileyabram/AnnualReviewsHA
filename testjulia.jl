# Chase Abram
# Used for testing how to run on Midway

task_id = string(ENV["SLURM_ARRAY_TASK_ID"])

io = open(task_id*"test.txt", "w");

write(io, "Test on "*task_id);

close(io);

println("Print statement for "*task_id*" again!")


# For Slurm
# module load julia

# julia testjulia.jl > /home/abram/AnnualReviewsHA/output