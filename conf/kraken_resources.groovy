/* 
Functions to calculate resources based on krakendb size
*/

// Function to calculate the folder size, excluding tar files
def calculateFolderSize(folderPath) {

    def folder = new File(folderPath)
    def totalSize = 0
    folder.eachFileRecurse { file ->
        if (file.isFile() && !file.name.endsWith('.tar.gz') && !file.name.endsWith('.tar')) {
            totalSize += file.length()
        }
    }
    return totalSize
}

// Function to calculate the memory required based on the folder size
def calculateMemoryForKraken(dbPath) {

    def folderSizeBytes = calculateFolderSize(dbPath)
    def folderSizeGB = folderSizeBytes / (1024 * 1024 * 1024)

    def memoryNeeded
    if (folderSizeGB < 1) {
        memoryNeeded = 2.GB
    } else {
        memoryNeeded = (folderSizeBytes / (1024 * 1024 * 1024)).toBigDecimal().GB + 5.GB 
    } // Adding 5 GB to the calculation for buffer

    log.info "Calculated memory needed: ${memoryNeeded}" // Printing out calculated memory
    return memoryNeeded
}