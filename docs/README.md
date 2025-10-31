# Krocus Documentation

## Overview

Krocus is a tool for rapid multi-locus sequence typing (MLST) directly from uncorrected long reads. This documentation provides comprehensive information for users, developers, and operators.

## Documentation Structure

### For Users

- **[USER_GUIDE.md](USER_GUIDE.md)** - Complete user manual
  - Installation instructions
  - Quick start guide
  - Detailed usage examples
  - Parameter explanations
  - Output interpretation
  - Troubleshooting

### For Operators

- **[WORK_INSTRUCTION.md](WORK_INSTRUCTION.md)** - Step-by-step procedures
  - Standard operating procedures
  - Quality control guidelines
  - Result interpretation criteria
  - Documentation requirements
  - Troubleshooting workflows

### For Developers

- **[API_DOCUMENTATION.md](API_DOCUMENTATION.md)** - Technical API reference
  - Module descriptions
  - Class and method documentation
  - Code examples
  - Testing guidelines
  
- **[DEVELOPER_GUIDE.md](DEVELOPER_GUIDE.md)** - Development guidelines
  - Architecture overview
  - Development setup
  - Coding standards
  - Testing procedures
  - Contribution guidelines

## Quick Links

### Getting Started

1. **Installation**: See [USER_GUIDE.md](USER_GUIDE.md#installation)
2. **Quick Start**: See [USER_GUIDE.md](USER_GUIDE.md#quick-start)
3. **First Analysis**: See [WORK_INSTRUCTION.md](WORK_INSTRUCTION.md#procedure)

### Common Tasks

- **Download database**: [USER_GUIDE.md](USER_GUIDE.md#krocus_database_downloader)
- **Run analysis**: [USER_GUIDE.md](USER_GUIDE.md#krocus)
- **Interpret results**: [WORK_INSTRUCTION.md](WORK_INSTRUCTION.md#step-5-interpret-results)
- **Real-time analysis**: [USER_GUIDE.md](USER_GUIDE.md#real-time-analysis)

### Reference

- **API Documentation**: [API_DOCUMENTATION.md](API_DOCUMENTATION.md)
- **Command reference**: [WORK_INSTRUCTION.md](WORK_INSTRUCTION.md#appendix-a-command-reference)
- **Troubleshooting**: [USER_GUIDE.md](USER_GUIDE.md#troubleshooting)

## Key Concepts

### Multi-Locus Sequence Typing (MLST)

MLST is a technique for characterizing bacterial isolates based on sequences of internal fragments of housekeeping genes. Each unique combination of alleles is assigned a sequence type (ST).

### K-mers

K-mers are subsequences of length k. Krocus uses k-mers to rapidly identify alleles without requiring read alignment or assembly. The k-mer size should match the error profile of your reads.

### Sequence Type (ST)

A sequence type is a unique identifier assigned to a specific combination of alleles across all MLST genes for a species.

### Coverage

Coverage refers to the percentage of expected k-mers that were found in the reads. Higher coverage indicates more confident results.

## Typical Workflow

```
1. Download Database
   ↓
2. Prepare FASTQ File
   ↓
3. Run Krocus
   ↓
4. Interpret Results
   ↓
5. Document Findings
```

## Performance Expectations

- **Speed**: 90 seconds average for typical bacterial genome
- **Memory**: ~1MB RAM per 1MB FASTQ input
- **Sensitivity**: 94% (in validation study)
- **Specificity**: 97% (in validation study)

## Supported Platforms

- **Sequencing**:
  - PacBio (RSII, Sequel, Sequel II)
  - Oxford Nanopore (MinION, GridION, PromethION)
  
- **Operating Systems**:
  - Linux (all distributions)
  - macOS
  - Windows (via WSL)

## System Requirements

### Minimum

- Python 3.3+
- 1GB RAM
- 100MB disk space (plus database size)

### Recommended

- Python 3.8+
- 4GB RAM
- 1GB disk space (plus database size)
- Multi-core processor for faster analysis

## Citation

If you use Krocus in your research, please cite:

**Krocus:**
> Andrew J. Page, Jacqueline A. Keane. (2018) Rapid multi-locus sequence typing direct from uncorrected long reads using Krocus. PeerJ 6:e5233 https://doi.org/10.7717/peerj.5233

**PubMLST:**
> Keith A. Jolley, James E. Bray, Martin C. J. Maiden. (2018) Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Research. 3:124 https://doi.org/10.12688/wellcomeopenres.14826.1

## Support

- **GitHub Issues**: https://github.com/andrewjpage/krocus/issues
- **Email**: andrewjpage+krocus@gmail.com
- **Documentation**: https://github.com/andrewjpage/krocus/tree/master/docs

## License

Krocus is licensed under GNU General Public License version 3 (GPL-3.0).

See the LICENSE file in the root directory for full license text.

## Contributing

Contributions are welcome! Please see [DEVELOPER_GUIDE.md](DEVELOPER_GUIDE.md) for:
- Code style guidelines
- Testing requirements
- Pull request process
- Development setup

## Version Information

- **Current Version**: See VERSION file in root directory
- **Changelog**: See CHANGELOG in root directory
- **Release Notes**: See GitHub releases page

## Additional Resources

### External Links

- **GitHub Repository**: https://github.com/andrewjpage/krocus
- **PyPI Package**: https://pypi.org/project/krocus/
- **PubMLST**: https://pubmlst.org/
- **Paper**: https://doi.org/10.7717/peerj.5233

### Related Tools

- **ARIBA**: Antimicrobial resistance identification from short reads
- **mlst**: Traditional MLST from assemblies
- **PubMLST BIGSdb**: Web-based MLST analysis

## Frequently Asked Questions

### Q: What read length is required?

A: Minimum 100bp, but longer reads (>1000bp) work better as they're more likely to span entire genes.

### Q: Can I use short reads (Illumina)?

A: Krocus is designed for long reads. For short reads, use traditional MLST tools like mlst or ARIBA.

### Q: How much coverage do I need?

A: Minimum 10x genome coverage recommended, 20x+ for best results.

### Q: Can I analyze mixed samples?

A: No, Krocus assumes a pure culture. Mixed samples will give unpredictable results.

### Q: How do I choose k-mer size?

A: Base it on your error rate. See [USER_GUIDE.md](USER_GUIDE.md#understanding-k-mer-size) for guidelines.

### Q: What if my species isn't available?

A: Check PubMLST.org - if a scheme exists there, you can download it with krocus_database_downloader.

### Q: Can I use my own database?

A: Yes, as long as it follows the PubMLST format (profile.txt and .tfa files).

## Updates and Maintenance

This documentation is maintained alongside the Krocus codebase. For the latest version, see the GitHub repository.

Last updated: Check git commit history for docs/ directory.
