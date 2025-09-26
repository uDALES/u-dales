# uDALES Variable Analysis - Report Generator

The `generate_reports.py` script generates comprehensive variable analysis reports for uDALES.

## Features

This script replaces the previous three separate scripts:
- ~~`generate_all_reports.py`~~ (moved to attic)
- ~~`generate_reports.py`~~ (moved to attic) 
- ~~`generate_status_report.py`~~ (moved to attic)

### Key Improvements:
- **Efficient processing**: Parses Fortran files only once, reuses data in memory for all reports
- **File tracking**: Full report now shows which file each namelist is defined in
- **Simplified interface**: No command-line arguments needed - automatically generates all reports

## Usage

```bash
# Generate both reports (no arguments needed)
python3 generate_reports.py
```

## Generated Reports

1. **`variable_status_report.md`** - 🎯 **Recommended starting point**
   - Color-coded overview with 🟢🟠🔴 status indicators
   - Executive summary with statistics
   - Categorized variables by namelist
   - Specific recommendations for improvements

2. **`variable_index_full.md`** - Comprehensive detailed report
   - Complete tables for each namelist
   - All variables with their status across 4 categories
   - Full variable listings and cross-references

## Quick Start

For a quick overview of variable status:
```bash
python3 generate_reports.py --report status
```

Then open `variable_status_report.md` to see the color-coded analysis.