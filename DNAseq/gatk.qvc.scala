 import java.io.File
import scala.util.Random

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.function.ListWriterFunction

class VariantDiscovery extends QScript {
  qscript =>
  
  /****************************************************************************
   * Required Parameters
   *****************************************************************************/
  
  @Input(doc="The reference file for the bam files.", shortName="R", fullName="reference_sequence")
  var referenceFile: File = _ // _ is scala shorthand for null
  
  @Input(doc="Bam file to relalign", shortName="I", fullName="input_file")
  var gvcfList: File = _

  @Output(doc = "Name of output", fullName="out")
  var out: File = _
  
  /****************************************************************************
   * Optional Parameters
   *****************************************************************************/
  
  @Argument(doc="scatter parameter", shortName="P", fullName="scatter_parameter", required=false)
  var scatter:  Int = _  

  @Argument(doc="nt parameter", shortName="N", fullName="num_threads", required=false)
  var nt:  Int = _
  
  @Argument(doc="nct parameter", shortName="C", fullName="num_cpu_threads_per_data_thread", required=false)
  var nct:  Int = _

  @Input(doc="Intervals to realign", shortName="L", required = false)
  var intervals: File = _

  /****************************************************************************
   * CommonArguments
   *****************************************************************************/
  
  trait CommonArguments extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFile
    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    this.memoryLimit = 20
  }
  
  
  /****************************************************************************               
   * Main script
   *****************************************************************************/
  
  def script() {
    
    val jointGVCF = new GenotypeGVCFs with CommonArguments
    val vqsrSNP = new VariantRecalibrator with CommonArguments
    val applySNPRecalibration = new ApplyRecalibration with CommonArguments
    val vqsrINDEL = new VariantRecalibrator with CommonArguments
    val applyINDELRecalibration = new ApplyRecalibration with CommonArguments

    jointGVCF.variant +:= qscript.gvcfList
    jointGVCF.out = qscript.out
    jointGVCF.dbsnp = new File("/scratch/ulg/genan/rgularte/mus_musculus/Ensembl/Grcm38/Annotation/Variation/Mus_musculus.sorted.vcf")
    jointGVCF.annotation = Seq("QualByDepth", "FisherStrand", "HaplotypeScore", "MappingQualityRankSumTest" , "ReadPosRankSumTest" , "VariantType")
    jointGVCF.nt = qscript.nt

    vqsrSNP.input +:= jointGVCF.out
    vqsrSNP.resource +:= new TaggedFile("/scratch/genan/rgularte/mus_musculus/Ensemb/Grcm38/Annotation/Variation/JAX-MDA1.vcf", "knwon=false,training=true,truth=true,prior=15.0")
    vqsrSNP.resource +:= new TaggedFile("/scratch/genan/rgularte/mus_musculus/Ensemb/Grcm38/Annotation/Variation/Mus_musculus_SNP_MO.vcf", "known=false,training=true,truth=false,prior=8.0")
    vqsrSNP.resource +:= new TaggedFile("/scratch/genan/rgularte/mus_musculus/Ensemb/Grcm38/Annotation/Variation/mgp.v3.snps.rsIDdbSNPv137.vcf.gz", "known=true,training=false,truth=false,prior=5.0")
    vqsrSNP.use_annotation = Seq("QD", "FS", "DP", "MQRankSum", "ReadPosRankSum")
    vqsrSNP.mode = org.broadinstitute.gatk.tools.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
    vqsrSNP.TStranche = Seq("99", "95", "90")
    vqsrSNP.recalFile = swapExt(jointGVCF.out, "Variants.vcf.gz",  "SNP.recal")
    vqsrSNP.tranchesFile = swapExt(jointGVCF.out, "Variants.vcf.gz", "SNP.tranches")
    vqsrSNP.rscriptFile = swapExt(jointGVCF.out, "Variants.vcf.gz", "SNP.plots.R")
    vqsrSNP.nt = qscript.nt

    applySNPRecalibration.input +:= jointGVCF.out
    applySNPRecalibration.out = swapExt(jointGVCF.out, "Variants.vcf.gz",  "SNP.vcf.gz")
    applySNPRecalibration.ts_filter_level = 99
    applySNPRecalibration.tranchesFile = vqsrSNP.tranchesFile
    applySNPRecalibration.recalFile = vqsrSNP.recalFile
    applySNPRecalibration.mode = org.broadinstitute.gatk.tools.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP
    applySNPRecalibration.nt = qscript.nt

    vqsrINDEL.input +:= jointGVCF.out
    vqsrINDEL.resource +:= new TaggedFile("/scratch/genan/rgularte/mus_musculus/Ensemb/Grcm38/Annotation/Variation/Mus_musculus_InDel_MO.vcf.gz", "knwon=false,training=true,truth=true,prior=10.0")
    vqsrINDEL.resource +:= new TaggedFile("/scratch/genan/rgularte/mus_musculus/Ensemb/Grcm38/Annotation/Variation/Mus_musculus.sorted.vcf.gz", "known=false,training=true,truth=false,prior=7.0")
    vqsrINDEL.resource +:= new TaggedFile("/scratch/genan/rgularte/mus_musculus/Ensemb/Grcm38/Annotation/Variation/mgp.v3.indels.rsIDdbSNPv137.vcf.gz", "known=true,training=false,truth=false,prior=5.0")
    vqsrINDEL.use_annotation = Seq("QD", "FS", "DP", "MQRankSum", "ReadPosRankSum")
    vqsrINDEL.mode = org.broadinstitute.gatk.tools.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
    vqsrINDEL.TStranche = Seq("99", "95", "90")
    vqsrINDEL.recalFile = swapExt(jointGVCF.out, "Variants.vcf.gz",  "INDEL.recal")
    vqsrINDEL.tranchesFile = swapExt(jointGVCF.out, "Variants.vcf.gz", "INDEL.tranches")
    vqsrINDEL.rscriptFile = swapExt(jointGVCF.out, "Variants.vcf.gz", "INDEL.plots.R")
    vqsrINDEL.nt = qscript.nt

    applyINDELRecalibration.input +:= jointGVCF.out
    applyINDELRecalibration.out = swapExt(jointGVCF.out, "Variants.vcf.gz",  "INDEL.vcf.gz")
    applyINDELRecalibration.ts_filter_level = 99
    applyINDELRecalibration.tranchesFile = vqsrINDEL.tranchesFile
    applyINDELRecalibration.recalFile = vqsrINDEL.recalFile
    applyINDELRecalibration.mode = org.broadinstitute.gatk.tools.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL
    applyINDELRecalibration.nt = qscript.nt


    add(jointGVCF, vqsrSNP, applySNPRecalibration, vqsrINDEL, applyINDELRecalibration)
    
  }
  
}	
