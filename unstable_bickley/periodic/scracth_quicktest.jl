
tmp = c_timeseries[100]
u1 = assemble(tmp).data[:, :, 1]
u2 = assemble(tmp).data[:, :, 2]
norm(u1-u2)