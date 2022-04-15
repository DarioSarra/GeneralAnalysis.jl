using Revise, GeneralAnalysis

data_dir = "/Volumes/GoogleDrive/My Drive/Flipping/Datasets/Stimulations/DRN_Opto_again"
bouts = CSV.read(joinpath(data_dir,"boutsDRN_Opto_again.csv"), DataFrame)
