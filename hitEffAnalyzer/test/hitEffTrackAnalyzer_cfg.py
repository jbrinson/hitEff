import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Demo")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

### Standard Configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration/Geometry/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
### For associating reco'd tracks to sim tracks
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")  
### Conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'DESIGN53_V15::All'

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

###Choose input file
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/w/wulsin/workspace/public/disappTrk/simStudyV2/CMSSW_5_3_3/src/muPartGun/batch/outputFEVT/SingleMuPt100_cfi_py_Ideal_RECO_FEVT_20.root'  
        )  
    )  

process.hitEff = cms.EDAnalyzer("hitEffAnalyzer",
    tracksTag    = cms.InputTag("generalTracks"), 
    tpTag        = cms.InputTag("mergedtruth", "MergedTrackTruth"),
    simtracksTag = cms.InputTag("g4SimHits"),
    )

###Output file name
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hitEff.root')
)

process.p = cms.Path(process.TrackRefitter + process.hitEff)

outfile = open('dumpedConfig.py','w'); print >> outfile,process.dumpPython(); outfile.close()

