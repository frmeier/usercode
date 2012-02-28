# Auto generated configuration file
# using: 
# Revision: 1.303.2.7 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: PYTHIA6_B0ToJpsiKs_7TeV_EvtGen_Frank_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW.root --step HLT:BPH,RAW2DIGI,L1Reco,RECO --beamspot Realistic7TeV2011Collision --conditions START42_V14A::All --pileup NoPileUp --datamix NODATAMIXER --eventcontent RECOSIM --datatier GEN-SIM-RECO --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('HLTrigger.Configuration.HLT_quarkonium_1E33_3E33_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('PYTHIA6_B0ToJpsiKs_7TeV_EvtGen_Frank_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.303.2.7 $'),
    annotation = cms.untracked.string('PYTHIA6_B0ToJpsiKs_7TeV_EvtGen_Frank_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW.root nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('PYTHIA6_B0ToJpsiKs_7TeV_EvtGen_Frank_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_root_HLT_quarkonium_RAW2DIGI_L1Reco_RECO.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)

# HLT stuff from https://espace.cern.ch/cms-quarkonia/onia-polarization/Generating%20Events/Detector%20simulation%20for%201E33.aspx
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1MuGMTParametersRcd' ),
        tag     = cms.string( 'L1MuGMTParameters_synctf_10_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
    ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1MuDTTFParametersRcd' ),
        tag     = cms.string( 'L1MuDTTFParameters_dttf11_TSC_09_17_col_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1MuCSCTFConfigurationRcd' ),
        tag     = cms.string( 'L1MuCSCTFConfiguration_90511_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1RPCBxOrConfigRcd' ),
        tag     = cms.string( 'L1RPCBxOrConfig_LHC7_1EX_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1RPCConeDefinitionRcd' ),
        tag     = cms.string( 'L1RPCConeDefinition_LHC7_1EX_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1RPCConfigRcd' ),
        tag     = cms.string( 'L1RPCConfig_LHC7_1EX_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1RPCHsbConfigRcd' ),
        tag     = cms.string( 'L1RPCHsbConfig_LHC7_1EX_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        )
    )

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START42_V14A::All'

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule()
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step])

