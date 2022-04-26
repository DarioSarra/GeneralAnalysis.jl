using Revise, GeneralAnalysis
##
data_dir = "/Volumes/GoogleDrive/My Drive/Flipping/Datasets/Stimulations/DRN_Opto_again"
streaks = CSV.read(joinpath(data_dir,"streaksDRN_Opto_again.csv"), DataFrame)
pokes = CSV.read(joinpath(data_dir,"pokesDRN_Opto_again.csv"), DataFrame)
bouts = CSV.read(joinpath(data_dir,"boutsDRN_Opto_again.csv"), DataFrame)
transform!(bouts,
    :Leave => (x-> contains.(x, "true")) => :Leave)
##
streaks.Pre_Interpoke = [contains(x,"NA") ? missing : parse(Float64,x) for x in streaks.Pre_Interpoke]
filter!(r -> r.Stim_Day &&
            r.Protocol_Day > 2 &&
            r.Streak_within_Block >2 &&
            !ismissing(r.Pre_Interpoke) &&
            r.Pre_Interpoke <=1,
            streaks)
##
streaks.Gen = levels!(categorical(streaks.Gen),["WT", "HET"])
isordered(streaks.Gen)
streaks.Protocol = levels!(categorical(streaks.Protocol),["90/90", "90/30","30/30"])
ordered!(streaks.Protocol, true)
streaks[9378,:Protocol] > streaks[1,:Protocol]
isordered(streaks.Protocol)
streaks.Stim = levels!(categorical(streaks.Stim),[false, true])
ordered!(streaks.Stim, true)
streaks.Wall = levels!(categorical(streaks.Wall),[false, true])
ordered!(streaks.Wall, true)
##
form = @formula(Leave ~ 1 + Omissions_plus_one + Stim&Gen + Protocol + Wall +
    (1+Protocol+Wall|MouseID));
form2 = @formula(Leave ~ 1 + Omissions_plus_one + (1|MouseID));
gm1 = GeneralizedLinearMixedModel(form, bouts, Bernoulli())
##
dyestuff = MixedModels.dataset(:dyestuff)
DataFrame(dyestuff)
##
fm = @formula(yield ~ 1 + (1|batch))
fm1 = fit(MixedModel, fm, dyestuff)
##
println.(propertynames(streaks))

##
df1 = combine(groupby(streaks,[:Stim,:Gen,:Protocol,:Wall]), :Num_pokes .=> [mean,sem])
@df df1 scatter(:Protocol.*:Gen.*string.(:Wall),:Num_pokes_mean, yerror = :Num_pokes_sem, group = :Stim,
    xrotation = 25)
##
form0 = @formula(AfterLast ~ 1  + Protocol + Wall +
    (1+Protocol+Wall|MouseID));
gm0 = fit(MixedModel,form0,streaks)
##
form1 = @formula(AfterLast ~ 1  + Protocol * Wall +
    (1+Protocol+Wall|MouseID));
gm1 = fit(MixedModel,form1,streaks)
mm_effectsize(gm0,gm1)
##
form2 = @formula(AfterLast ~ 1  + Protocol * Wall + Stim&Gen +
    (1+Protocol+Wall|MouseID));
gm2 = fit(MixedModel,form2,streaks)
mm_effectsize(gm1,gm2)
##
form3 = @formula(AfterLast ~ 1  + Protocol * Wall * Stim&Gen +
    (1+Protocol+Wall|MouseID));
gm3 = fit(MixedModel,form3,streaks)
mm_effectsize(gm1,gm3)
##
form4 = @formula(AfterLast ~ 1  + Protocol * Wall + Stim&Gen + Wall&Stim&Gen + Protocol&Stim&Gen +
    (1+Protocol+Wall|MouseID));
gm4 = fit(MixedModel,form4,streaks)
mm_effectsize(gm1,gm4)
##
form5 = @formula(AfterLast ~ 1  + Protocol * Wall *Stim * Gen+
    (1+Protocol+Wall|MouseID));
gm5 = fit(MixedModel,form5,streaks)
mm_effectsize(gm1,gm5)
##
streaks.Pre_Interpoke = [contains(x,"NA") ? missing : parse(Float64,x) for x in streaks.Pre_Interpoke]
df = filter(r-> r.Gen =="HET" && r.Wall &&
    r.Streak_within_Block >2 &&
    !ismissing(r.Pre_Interpoke) &&
    r.Pre_Interpoke <=1, streaks)

df1 = combine(groupby(df,[:Stim,:MouseID]),:AfterLast => mean => :AfterLast)
df1[!,:Stimulation] = [x ? "On" : "Off" for x in df1.Stim]
# open_html_table(df1)
df3 = unstack(df1[:,Not(:Stim)],:Stimulation, :AfterLast)
dropmissing!(df3)
es = CohenD(df3.On, df3.Off, quantile=0.95)
effectsize(es)
es = GlassΔ(df3.On, df3.Off, quantile=0.95)
effectsize(es)
es = HedgeG(df3.On, df3.Off, quantile=0.95)
effectsize(es)

df2 = combine(groupby(df1,[:Stim]),:AfterLast .=> [mean,std,sem])
union(df.MouseID)

√((SD12 + SD22) ⁄ 2)
(df2[2,:AfterLast_mean] - df2[1,:AfterLast_mean])/(√(((df2[2,:AfterLast_std])^2 + (df2[1,:AfterLast_std])^2)/2))

##
function meads_resource(T,B; E = 20)
    # N = total number of individuals or units in the study (minus 1)
    # E = degrees of freedom of the error component, and should be somewhere between 10 and 20
    # B = Blocking component - 1
    # T = treament levels -1
    # E = N - B - T
    #E = N*(T) - T - B - 3
    #N = (E + T + B + 3)/ (T)
    (E + T + B + 3)/ (T)
end
meads_resource(2,12)
##
rng = MersenneTwister(1234321)
nsims = 1000
samp = parametricbootstrap(rng, 1000, gm4);
df0 = DataFrame(samp.allpars);
df1 = DataFrame(samp.coefpvalues);
ptbl = power_table(samp, 0.05)
pretty_table(ptbl)
## Simulate
subj_n = 20
item_n = 500
subj_btwn = Dict(:gen => ["wt", "het"])
item_btwn = nothing
both_win = Dict(:stim => ["false", "true"], :patch => ["r", "m", "p"], :travel => ["l", "s"])

fake_data = simdat_crossed(subj_n, item_n,
                           subj_btwn = subj_btwn,
                           item_btwn = item_btwn,
                           both_win = both_win);

fake_data_df = DataFrame(fake_data)
open_html_table(fake_data_df)
contrasts = Dict(:stim => HelmertCoding(),
                 :patch => HelmertCoding(),
                 :travel => HelmertCoding(),
                 :gen => HelmertCoding(),
                 :item => Grouping(),
                 :subj => Grouping());

formula = @formula(dv ~ 1  + patch * travel + stim&gen + travel&stim&gen + patch&stim&gen +
     (1+patch+travel|subj));
formula = @formula(dv ~ 1  + patch * travel * stim * gen + (1+patch+travel|subj))
m1 = fit(MixedModel, formula, fake_data_df)

parametricbootstrap(rng, 1000, m1,use_threads = false)
sim1 = parametricbootstrap(rng, nsims, m1;
                           β = gm4.β,
                           σ = gm4.θ,
                           θ = gm4.θ,
                           use_threads = false);

##
# item response
println(VarCorr(gm4).σρ)
VarCorr(gm4)
corr_exact = VarCorr(gm4).σρ.item.ρ[1]
σ_residuals_exact = gm4.σ
σ_1_exact = VarCorr(gm4).σρ.item.σ[1] / σ_residuals_exact
σ_2_exact = VarCorr(gm4).σρ.item.σ[2] / σ_residuals_exact

re_item_corr = [1.0 corr_exact; corr_exact 1.0]
re_item = create_re(σ_1_exact, σ_2_exact; corrmat = re_item_corr)
# subject response
σ_residuals_exact = gm4.σ
σ_3_exact = VarCorr(gm4).σρ.MouseID.σ[1] / σ_residuals_exact
re_subj = create_re(σ_3_exact)
#
x = 36
x*0.1 + x
x - x*0.1
36/(2*3*2*2)
##
wf_bouts = CSV.read("/Users/dariosarra/Documents/Lab/Walton/WaltonForaging/DAphotometry/Processed/AllBouts.csv", DataFrame)


println.(propertynames(wf_bouts))
wf_dscr = combine(groupby(wf_bouts,[:Startdate, :SubjectID]),
    :Patch => maximum => :Patch,
    :Rewarded => (x -> sum(skipmissing(x)))=> :Rewarded
    )
mean(wf_dscr.Rewarded)
std(wf_dscr.Rewarded)
open_html_table(wf_dscr)
##
characteristic = ["small", "big"]
treatment = [false, true]
condition = ["low", "medium", "high"]
n_itmes = 100
group_size = 10

ex1 = DataFrame()
for c in characteristic
    for g in 1:group_size
        for i in 1:n_itmes
            effect = c == "big" ? 2 : 1
            res = rand() + effect + g * rand()
            push!(ex1, (Characteristic = c, Subject = string(g),
                Item = i, Response = res))
        end
    end
end
ex2 = copy(ex1)
ex2.Characteristic = categorical(ex2.Characteristic)
f = @formula(Response ~ 1 + Characteristic + (1+Characteristic|Subject))
ex1_m = fit(MixedModel,f,ex1)
ex2_m = fit(MixedModel,f,ex2)
##
streaks = CSV.read(joinpath(data_dir,"streaksDRN_Opto_again.csv"), DataFrame)
streaks.Pre_Interpoke = [contains(x,"NA") ? missing : parse(Float64,x) for x in streaks.Pre_Interpoke]
filter!(r -> r.Stim_Day &&
            r.Protocol_Day > 2 &&
            r.Streak_within_Block >2 &&
            !ismissing(r.Pre_Interpoke) &&
            r.Pre_Interpoke <=1,
            streaks)
##
ex_v1 = streaks[:,[:AfterLast, :Protocol,:Wall,:Gen,:Stim,:MouseID]]
rename!(ex_v1, :AfterLast => :Response, :MouseID => :SubjectID)
CSV.write(joinpath("/Users/dariosarra/Desktop","Example.csv"), ex_v1)
describe(ex_v1)
ex_v2 = copy(ex_v1)
for x in [:Protocol,:Wall,:Gen,:Stim]
    ex_v2[!,x] = categorical(ex_v2[:,x])
end
levels!(ex_v2.Gen,["WT", "HET"])
ex_v3 = copy(ex_v2)
levels!(ex_v3.Protocol,["90/90", "90/30","30/30"])
isordered(ex_v3.Protocol)
contrasts = Dict(:Protocol => EffectsCoding(; base="90/90", levels = ["90/90", "90/30","30/30"]),
                        :Gen => EffectsCoding(; base="WT"),
                        :Stim => DummyCoding(; base=false),
                        :Wall => DummyCoding(; base=false),
                        :SubjectID => Grouping())
contrasts = Dict(:Protocol => EffectsCoding(; base="90/90"),
                        :Gen => EffectsCoding(; base="WT"),
                        :Stim => EffectsCoding(; base=false),
                        :Wall => EffectsCoding(; base=false),
                        :SubjectID => Grouping())
form = @formula(Response ~ 1  + Protocol * Wall *Stim * Gen+
    (1+Protocol+Wall|SubjectID));
m_v1 = fit(MixedModel,form,ex_v1;contrasts)
m_v2 = fit(MixedModel,form,ex_v2;contrasts)
m_v3 = fit(MixedModel,form,ex_v3;contrasts)
