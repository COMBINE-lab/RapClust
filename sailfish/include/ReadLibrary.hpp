#ifndef READ_LIBRARY_HPP
#define READ_LIBRARY_HPP

#include <vector>
#include <exception>
#include <set>

#include <boost/filesystem.hpp>

#include "LibraryFormat.hpp"

/**
 * This class represents the basic information about a library of reads, like
 * its paired-end status, the reads that should appear on the forward and reverse strand,
 * and the relative orientation of the reads.
 */
class ReadLibrary {
public:
    /**
     * Construct a new ReadLibrary of the given format
     */
    ReadLibrary(LibraryFormat& fmt) : fmt_(fmt),
        libTypeCounts_(std::vector<std::atomic<uint64_t>>(LibraryFormat::maxLibTypeID() + 1)){}

    /**
     * Copy constructor
     */
    ReadLibrary(const ReadLibrary& rl) :
        fmt_(rl.fmt_),
        unmatedFilenames_(rl.unmatedFilenames_),
        mateOneFilenames_(rl.mateOneFilenames_),
        mateTwoFilenames_(rl.mateTwoFilenames_),
        libTypeCounts_(std::vector<std::atomic<uint64_t>>(LibraryFormat::maxLibTypeID() + 1)) {
            auto mc = LibraryFormat::maxLibTypeID() + 1;
            for (size_t i = 0; i < mc; ++i) { libTypeCounts_[i].store(rl.libTypeCounts_[i].load()); }
        }

    /**
     * Move constructor
     */
    ReadLibrary(ReadLibrary&& rl) :
        fmt_(rl.fmt_),
        unmatedFilenames_(std::move(rl.unmatedFilenames_)),
        mateOneFilenames_(std::move(rl.mateOneFilenames_)),
        mateTwoFilenames_(std::move(rl.mateTwoFilenames_)),
        libTypeCounts_(std::vector<std::atomic<uint64_t>>(LibraryFormat::maxLibTypeID() + 1)) {
            auto mc = LibraryFormat::maxLibTypeID() + 1;
            for (size_t i = 0; i < mc; ++i) { libTypeCounts_[i].store(rl.libTypeCounts_[i].load()); }
        }

    /**
     * Add files containing mated reads (from pair 1 of the mates) to this library.
     */
    void addMates1(const std::vector<std::string>& mateOneFilenames) {
        mateOneFilenames_ = mateOneFilenames;
    }

    /**
     * Add files containing mated reads (from pair 2 of the mates) to this library.
     */
    void addMates2(const std::vector<std::string>& mateTwoFilenames) {
        mateTwoFilenames_ = mateTwoFilenames;
    }

    /**
     * Add files containing unmated reads.
     */
    void addUnmated(const std::vector<std::string>& unmatedFilenames) {
        unmatedFilenames_ = unmatedFilenames;
    }

    /**
     * Return true if this read library is for paired-end reads and false otherwise.
     */
    bool isPairedEnd() {
        return (fmt_.type == ReadType::PAIRED_END);
    }


    bool checkFileExtensions_(std::vector<std::string>& filenames, std::stringstream& errorStream) {
        namespace bfs = boost::filesystem;

        std::set<std::string> acceptableExensions = {".FASTA", ".FASTQ", ".FA", ".FQ",
                                                     ".fasta", ".fastq", ".fa", ".fq"};
        bool extensionsOK{true};
        for (auto& fn : filenames) {
            auto fpath = bfs::path(fn);
            auto ext = fpath.extension().string();
            if (bfs::is_regular_file(fpath)) {
                if (acceptableExensions.find(ext) == acceptableExensions.end()) {
                    errorStream << "ERROR: file [" << fn << "] has extension " << ext << ", "
                        << "which suggests it is neither a fasta nor a fastq file.\n"
                        << "Is this a compressed file?  If so, consider replacing: \n\n"
                        << fn << "\n\nwith\n\n"
                        << "<(decompressor " << fn << ")\n\n"
                        << "which will decompress the reads \"on-the-fly\"\n\n";
                    extensionsOK = false;
                } else if (bfs::is_empty(fpath)) {
                    errorStream << "ERROR: file [" << fn << "] appears to be empty "
                        "(i.e. it has size 0).  This is likely an error. "
                        "Please re-run salmon with a corrected input file.\n\n";
                        extensionsOK = false;
                }
            }
        }
        return extensionsOK;
    }

    bool isRegularFile() {
        if (isPairedEnd()) {
            for (auto& m1 : mateOneFilenames_) {
                if (!boost::filesystem::is_regular_file(m1)) { return false; }
            }
            for (auto& m2: mateTwoFilenames_) {
                if (!boost::filesystem::is_regular_file(m2)) { return false; }
            }
        } else {
            for (auto& um : unmatedFilenames_) {
                if (!boost::filesystem::is_regular_file(um)) { return false; }
            }
        }
        return true;
    }

    std::string readFilesAsString() {
        std::stringstream sstr;
        if (isPairedEnd()) {
            size_t n1 = mateOneFilenames_.size();
            size_t n2 = mateTwoFilenames_.size();
            if (n1 == 0 or n2 == 0 or n1 != n2) {
                sstr << "LIBRARY INVALID --- You must provide #1 and #2 mated read files with a paired-end library type";
            } else {
                for (size_t i = 0; i < n1; ++i) {
                    sstr << "( " << mateOneFilenames_[i] << ", " <<
                            mateTwoFilenames_[i] << " )";
                    if (i != n1 - 1) { sstr << ", "; }
                }
            }
        } else { // single end
            size_t n = unmatedFilenames_.size();
            if (n == 0) {
                sstr << "LIBRARY INVALID --- You must provide unmated read files with a single-end library type";
            } else {
                for (size_t i = 0; i < n; ++i) {
                    sstr << unmatedFilenames_[i];
                    if (i != n - 1) { sstr << ", "; }
                }
            }
        } // end else
        return sstr.str();
    }

    /**
     * Checks if this read library is valid --- if it's paired-end, it should have mate1/2 reads and the same
     * number of files for each; if it's unpaired it should have only unpaired files.
     * NOTE: This function throws an exception if this is not a valid read library!
     */
    void checkValid() {

        bool readsOK{true};

        std::stringstream errorStream;
        errorStream << "\nThe following errors were detected with the read files\n";
        errorStream << "======================================================\n";

        if (isPairedEnd()) {
            size_t n1 = mateOneFilenames_.size();
            size_t n2 = mateTwoFilenames_.size();
            if (n1 == 0 or n2 == 0 or n1 != n2) {
                errorStream << "You must provide #1 and #2 mated read files with a paired-end library type\n";
                readsOK = false;
            }
        } else {
            size_t n = unmatedFilenames_.size();
            if (n == 0) {
                errorStream << "You must provide unmated read files with a single-end library type\n";
                readsOK = false;
            }

        }

        // Check if the user tried to pass in non-fast{a,q} files.  If so,
        // throw an exception with the appropriate error messages.
        // NOTE: This check currently does nothing useful for non-regular files
        // (i.e. named-pipes).  If the user passed in a non-regular file, we should
        // have some other mechanism to check if it's of an expected format and provide
        // a reasonable error message otherwise.
        readsOK = readsOK && checkFileExtensions_(mateOneFilenames_, errorStream);
        readsOK = readsOK && checkFileExtensions_(mateTwoFilenames_, errorStream);
        readsOK = readsOK && checkFileExtensions_(unmatedFilenames_, errorStream);

        if (!readsOK) {
            throw std::invalid_argument(errorStream.str());
        }
    }

    /**
     * Return the vector of files containing the mate1 reads for this library.
     */
    const std::vector<std::string>& mates1() const { return mateOneFilenames_; }

    /**
     * Return the vector of files containing the mate2 reads for this library.
     */
    const std::vector<std::string>& mates2() const { return mateTwoFilenames_; }

    /**
     * Return the vector of files containing the unmated reads for the library.
     */
    const std::vector<std::string>& unmated() const { return unmatedFilenames_; }

    /**
     * Return the LibraryFormat object describing the format of this read library.
     */
    const LibraryFormat& format() const { return fmt_; }

    /**
    * Update the library type counts for this read library given the counts
    * in the vector `counts` which has been passed in.
    */
    inline void updateLibTypeCounts(const std::vector<uint64_t>& counts) {
        size_t lc{counts.size()};
        for (size_t i = 0; i < lc; ++i) { libTypeCounts_[i] += counts[i]; }
    }

    std::vector<std::atomic<uint64_t>>& libTypeCounts() {
        return libTypeCounts_;
    }

private:
    LibraryFormat fmt_;
    std::vector<std::string> unmatedFilenames_;
    std::vector<std::string> mateOneFilenames_;
    std::vector<std::string> mateTwoFilenames_;
    std::vector<std::atomic<uint64_t>> libTypeCounts_;
};

#endif // READ_LIBRARY_HPP
