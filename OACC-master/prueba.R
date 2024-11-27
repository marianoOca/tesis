
#Directories are changed to work from phython
source("scripts/BDM1D.R")
source("scripts/compressionLength.R")

my_function <- function(bdmInputString) {
  bdmAlphabet <- 256
  blockSize <- 12
  blockOverlap <- 0

  values <- c ()

  if (bdmAlphabet == 256) {

    # convert UTF-8 string to binary
    binString <- getBinString(bdmInputString)

    values[1] <- paste0(
      sprintf("%.4f",stringBDM(
        splitString(binString,
                    blockSize = blockSize, 
                    offset = blockSize - blockOverlap),
                    base = 2)),
      " bits")

    values[2] <- paste0(
      sprintf("%.4f",stringBDMLD(
        splitString(binString,
                    blockSize = blockSize, 
                    offset = (blockSize - blockOverlap)),
                    base = bdmAlphabet)), 
      " steps")

  } else {

    values[1] <- paste0(
      sprintf("%.4f",stringBDM(
        splitString(bdmInputString,
                    blockSize = blockSize, 
                    offset = blockSize - blockOverlap),
                    base = bdmAlphabet)), 
      " bits")
  }

  if (bdmAlphabet == 2){
    values[2] <- paste0(
      sprintf("%.4f",stringBDMLD(
        splitString(bdmInputString,
                    blockSize = blockSize,
                    offset = blockSize -blockOverlap),
                    base = bdmAlphabet)),
      " steps")
  }

  #entropy
  values[3] <- paste0(sprintf("%.4f",
                              entropy(bdmInputString)),
                      " bit(s)")

  #second order entropy
  values[4] <- paste0(sprintf("%.4f",
                              entropy2(bdmInputString)),
                      " bit(s)")

  #compression length
  values[5] <- paste0(compressionLength(bdmInputString,
                                        "gzip") * 8,
                      " bits")

  values[6] <- nchar(bdmInputString)
  values[7] <- countSymbols(bdmInputString)
  values[8] <- bdmAlphabet
  values[9] <- blockSize
  values[10] <- blockOverlap

  resultRowNames  <-  c("BDM algorithmic complexity estimation",
                        "BDM logical depth estimation",
                        "Shannon entropy",
                        "Second order entropy",
                        "Compression length (using gzip)",
                        "String length",
                        "# of symbols in string",
                        "# of symbols in CTM alphabet",
                        "Block size",
                        "Block overlap")

  if (!(bdmAlphabet == 2 || bdmAlphabet == 256)) {
    values <- values[-2]
    resultRowNames <- resultRowNames[-2]
  }

  result <- data.frame(values)
  rownames(result) <- resultRowNames

  result
}

my_function("MAKFVKVKGTDQVLVKKEISILNIARHRNILHLHESFESMEELVMIFEFISGLDIFERINTSAFELNEREIVSYVHQVCEALQFLHSHNIGHFDIRPENIIYQTRRSSTIKIIEFGQARQLKPGDNFRLLFTAPEYYAPEVHQHDVVSTATDMWSLGTLVYVLLSGINPFLAETNQQIIENIMNAEYTFDEEAFKEISIEAMDFVDRLLVKERKSRMTASEALQHPWLKQKIERVSTKVIRTLKHRRYYHTLIKKDLNMVVSAARISCGGAIRSQKGVSVAKVKVASIEIGPVSGQIMHAVGEEGGHVKYVCKIENYDQSTQVTWYFGVRQLENSEKYEITYEDGVAILYVKDITKLDDGTYRCKVVNDYGEDSSYAELFVKGVREVYDYYCRRTMKKIKRRT")