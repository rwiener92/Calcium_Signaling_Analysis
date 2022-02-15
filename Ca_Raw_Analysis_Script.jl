# This Julia code is designed to take Results (as CSV data) from FIJI analysis of Ca2+ Flux experiments
# This code will process raw CSV data and return
#       1. Raw signal for inspection and filtering
#       2. Raw signal normalized to baseline (see ΔS formula)
#       3. An average signal for all normalized cells
#       4. A replotted subplot of any filtered cells
# Inputs: cd, Tx_cut, remove_cell, base_stop, CSV.write
### Written by Robert J. Wiener (c) 09/28/2021
##
using JLD
using CSV
using DataFrames
using Statistics
using Plots; plotlyjs()

# Change directory to exp folder containing Results.csv from FIJI output
#cd("C:\\Users\\Rob\\Documents\\Costa Lab\\RESEARCH PROJECTS\\MFS hiPSC-ASMCs (RJW PhD)\\Calcium Signaling");
pwd()
#########################################################################

#### DATA IMPORT ####
#####################
# CSV import FIJI Results as DataFrame, then convert DataFrame to Matrix
df = CSV.read("Results.csv", DataFrame);
matdat = Matrix(df);

#### CREATE MATRIX & SPECIFY VARS ####
######################################
# Isolate time data (col_1), Background (col_2), Cell Signal (col_3-end)
Time = matdat[:,1];
Background = matdat[:,2];
Cellsig = matdat[:,3:end];
Tx_cut = size(Cellsig, 1); #specify start time for TritonX additon (use 'size(Cellsig, 1)' if no Tx)


#### CONDITIONAL FILTERING ####
###############################
### Input cell number to remove ###
remove_cell = [:x2]; #keep vector empty for none, otherwise index with symbol e.g. :x8 (for cell 8)
if isempty(remove_cell) # Starting run i.e. no cells removed

### Plot Raw Data vs. Time
# inspect cell signals and confirm no abrupt changes from fluid ctrl (@30 sec)
fig_raw = plot(Time, Background, label="Background");
    for i in 1:size(Cellsig, 2)
        fig_raw = plot!(Time, Cellsig[:,i], label="Cell0$(i)", legend=:outertopright)
    end
    title!("Raw")
    xlabel!("Time (sec)")
    ylabel!("Raw Fluo.")
# View raw data and consider any cells to remove
display(fig_raw)
   
else # Moves here if cells are being removed

    # Convert Cellsig matrix back to DF & remove select cells, then return Cellsig to Matrix
    temprem = DataFrame(Cellsig, :auto);
    select!(temprem, Not(remove_cell));
    Cellsig = Matrix(temprem);
    #= depreciated loop removal (bug in OOR selection)
    for i in 1:length(remove_cell)
        Cellsig = Cellsig[:, 1:size(Cellsig, 2) .!= remove_cell[i]]
    end
    =#

    fig_raw_rem = plot(Time, Background, label="Background");
    for i in 1:size(Cellsig, 2)
        fig_raw_rem = plot!(Time, Cellsig[:,i], label="Cell0$(i)_rawrem", legend=:outertopright)
    end
    title!("Raw Rem")
    xlabel!("Time (sec)")
    ylabel!("Raw Fluo.")
# View new raw data with cells removed
display(fig_raw_rem)
end



#### NORMALIZATION ####
#######################
# Use this to cut fluid ctrl time (i.e. row 31-59)
#Cellsig = matdat[[1:30; 60:end], :]

### Calculate normalized change in signal per cell (ΔS_cell)
# use 1-60 for no cut fluid ctrl time, use 1-30 for cut
base_stop = 59; #sec
Base_Signal = mean(Background[1:base_stop]);

ΔS_cell = zeros(size(Cellsig))
for j in 1:size(Cellsig, 2), i in 1:size(Cellsig, 1)
ΔS_cell[i, j] = (Cellsig[i, j] - Background[i]) / (mean(Cellsig[1:base_stop, j]) - Base_Signal)
end


#### PLOTTING ####
##################
### Plot ΔS_cell vs. Time
fig_ΔS = plot();
for i in 1:size(ΔS_cell, 2)
fig_ΔS = plot!(Time, ΔS_cell[1:Tx_cut,i], label="ΔS_Cell0$i", legend=:outertopright)
end
title!("Normalized Signal")
xlabel!("Time (sec)")
ylabel!("RFU.")
display(fig_ΔS)

# Average Cell signal
ΔS_cell_mean = mean(ΔS_cell, dims=2);
fig_ΔS_cell_mean = plot(Time, ΔS_cell_mean[1:Tx_cut], label="ΔS_cell_mean [$(size(Cellsig, 2)) cells]", legend=:outertopright);
title!("Average Normalized Signal")
xlabel!("Time (sec)")
ylabel!("RFU.")
display(fig_ΔS_cell_mean)



#########################
### Final Subplotting ###
if isempty(remove_cell)
    plot(fig_raw, fig_ΔS, fig_ΔS_cell_mean)
else
    plot(fig_raw, fig_raw_rem, fig_ΔS, fig_ΔS_cell_mean)
end

#CSV.write("MFS-PM, S6.csv", DataFrame(ΔS_cell_mean, :auto))

