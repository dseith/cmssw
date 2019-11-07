seed=$(($(date +%s) % 10000 + 1))
nevents=10
echo $seed

source settings18.txt

echo $configuration, $fileout_gen, $fileout_premix, $fileout_aod, $fileout_miniaod, $conditions, $beamspot, $era, $python_filename_gen, $python_filename_premix, $pileup_input
echo "########################"

#########
#LHE, GEN
#########
cmsDriver.py $configuration --fileout $fileout_gen  --mc --eventcontent RAWSIM,LHE --datatier GEN,LHE --conditions $conditions  --beamspot $beamspot  --step LHE,GEN --geometry DB:Extended --era $era --python_filename $python_filename_gen  --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring --customise_commands process.RandomNumberGeneratorService.externalLHEProducer.initialSeed="int(${seed})" -n $nevents --nThreads 20 

# #########
# #SIM
# #########
# cmsDriver.py step1 --filein $fileout_gen --fileout $fileout_sim --mc --eventcontent RAWSIM --runUnscheduled --datatier SIM --conditions $conditions --beamspot $beamspot --step SIM --nThreads 8 --geometry DB:Extended --era $era --python_filename $python_filename_sim --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n $nevents
# 
# #########
# #PREMIX (DIGI, DATAMIX, L1, DIGI2RAW)
# #########
# cmsDriver.py step1 --filein $fileout_sim --fileout $fileout_premix  --pileup_input $pileup_input --mc --eventcontent PREMIXRAW --runUnscheduled --datatier GEN-SIM-DIGI --conditions $conditions --step DIGI,DATAMIX,L1,DIGI2RAW --procModifiers premix_stage2 --nThreads 8 --geometry DB:Extended --datamix PreMix --era $era --python_filename $python_filename_premix --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n $nevents
# 
# #########
# #HLT #conditions are somehow different??? Maybe one has to use the conditions of the triggers??
# #########
# cmsDriver.py step1 --filein $fileout_premix --fileout $fileout_hlt --mc --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions $conditions --customise_commands 'process.source.bypassVersionCheck = cms.untracked.bool(True)' --step $hltstep --nThreads 8 --geometry DB:Extended --era $era --python_filename $python_filename_hlt --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n $nevents
# 
# #########
# #AOD/RECO
# #########
# cmsDriver.py step1 --filein $fileout_hlt --fileout $fileout_aod --mc --eventcontent AODSIM --runUnscheduled --datatier AODSIM --conditions $conditions --step RAW2DIGI,L1Reco,RECO,RECOSIM --nThreads 16 --geometry DB:Extended --era $era --python_filename $python_filename_aod --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n $nevents
# 
# 
# #########
# #MINIAOD
# #########
# cmsDriver.py step1 --filein $fileout_aod --fileout $fileout_miniaod --mc --eventcontent MINIAODSIM --runUnscheduled --datatier MINIAODSIM --conditions $conditions --step PAT --nThreads 8 --geometry DB:Extended --era $era --python_filename $python_filename_miniaod --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n $nevents
# 
# 
# #########
# #NANOAOD
# #########
# cmsDriver.py step1 --filein $fileout_miniaod --fileout $fileout_nanoaod --mc --eventcontent NANOEDMAODSIM --datatier NANOAODSIM --conditions $conditions --step NANO --nThreads 8 --era $era --python_filename $python_filename_nanoaod --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n $nevents
# 
# 
# 
# 
