import java.io.File
import scala.util.Random

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.function.ListWriterFunction

// package org.broadinstitute.gatk.queue.qscripts

class HapCaller extends QScript {
  qscript =>
  
  /****************************************************************************
   * Required Parameters
   *****************************************************************************/
  
  @Input(doc="The reference file for the bam files.", shortName="R", fullName="reference_sequence")
  var referenceFile: File = _ // _ is scala shorthand for null
  
  @Input(doc="Bam file to relalign", shortName="I", fullName="input_file")
  var myBam: File = _
  
  
  /****************************************************************************
   * Optional Parameters
   *****************************************************************************/
  
  @Argument(doc="scatter parameter", shortName="P", fullName="scatter_parameter", required=false)
  var scatter:  Int = _  
  
  @Argument(doc="nt parameter", shortName="N", fullName="num_threads", required=false)
  var nt:  Int = _
  
  @Argument(doc="nct parameter", shortName="C", fullName="num_cpu_threads_per_data_thread", required=false)
  var nct:  Int = _
  
//  @Argument(doc="genotype mode", shortName="gt_mode", required = false)
//  var gt_mode:  String = _
  
//  @Argument(doc="genotype mode", shortName="pcrModel", fullName = "pcr_indel_model", required = false)
//  var pim:  String = _

  @Argument(doc="max alleles", shortName="maxAltAlleles", fullName = "max_alternate_alleles", required = false)
  var maxAltAlleles:  Int = _

//  @Argument(doc="ERC", shortName="ERC", required = false)
//  var ERC:  String = _

//  @Argument(doc="variant index type", fullName="variant_index_type", required = false)
//  var vit:  String = _

//  @Argument(doc="variant index parameter", fullName="variant_index_parameter", required = false)
//  var vip:  Int = _

//  @Argument(doc="pair hidden markov model", shortName="pairHMM", required = false)
//  var pairHMM:  String = _

  @Input(doc="dbSNP file", shortName="D", fullName="dbsnp", required=false)
  var dbsnp: File = _
  
  @Input(doc="Intervals to realign", shortName="L", required = false)
  var intervals: File = _
  
  
  /****************************************************************************
   * Classes (GATK Walkers)
   *****************************************************************************/
  
  trait HCArguments extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFile
    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    this.memoryLimit = 12
//    this.variant_index_parameter = 128000
  }
  
    /****************************************************************************               
   * Main script
   *****************************************************************************/
  
  def script() {
    
    val hc = new HaplotypeCaller with HCArguments
    
    hc.scatterCount = qscript.scatter
    hc.input_file +:= qscript.myBam
//    hc.out = swapExt(qscript.myBam, ".dd.ir.bqsr.bam", ".vcf")
    hc.nct = qscript.nct
    hc.pairHMM = org.broadinstitute.gatk.utils.pairhmm.PairHMM.HMM_IMPLEMENTATION.VECTOR_LOGLESS_CACHING
    hc.pcr_indel_model = org.broadinstitute.gatk.tools.walkers.haplotypecaller.PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL.CONSERVATIVE
    hc.out = swapExt(qscript.myBam, ".dd.ir.bqsr.bam", ".g.vcf")
    hc.emitRefConfidence = org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode.GVCF
    hc.variant_index_type = org.broadinstitute.gatk.utils.variant.GATKVCFIndexType.LINEAR
    hc.variant_index_parameter = 128000 

    add(hc)
    
  }
  
}
