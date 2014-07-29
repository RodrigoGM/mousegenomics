import java.io.File
import scala.util.Random

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.function.ListWriterFunction

class IndelRealignmet extends QScript {
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
  
  @Input(doc="dbSNP file", shortName="D", fullName="dbsnp", required=false)
  var dbsnp: File = _

  @Input(doc="Known in/del VCF file", shortName="known", required=false)
  var known: File = _

  @Input(doc="Intervals to realign", shortName="L", required = false)
  var intervals: File = _
  
  
  /****************************************************************************
   * Classes (GATK Walkers)
   *****************************************************************************/
  
  trait CommonArguments extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFile
    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    this.memoryLimit = 12
  }
  
  
  /****************************************************************************               
   * Main script
   *****************************************************************************/
  
  def script() {
    
      val targetCreator = new RealignerTargetCreator with CommonArguments
      val indelRealigner = new IndelRealigner with CommonArguments
      
      targetCreator.input_file +:= qscript.myBam
      targetCreator.out = swapExt(myBam, "bam", "intervals")
      targetCreator.nt = qscript.nt
      targetCreator.known +:= qscript.known

      indelRealigner.input_file +:= qscript.myBam
      indelRealigner.targetIntervals = targetCreator.out
      indelRealigner.out = swapExt(myBam, "dd.bam", "dd.ir.bam")
      indelRealigner.known +:= qscript.known

      add(targetCreator, indelRealigner)
      
  }
  
}	
