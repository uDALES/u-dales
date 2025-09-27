# uDALES Variable Analysis - Report Generator

The `test_sourcecode.py` script generates comprehensive variable analysis reports for uDALES (renamed from `generate_reports.py`).

## Usage

```bash
# Generate both reports (no arguments needed)
python3 test_sourcecode.py
```

## Generated Reports

1. **`test_sourcecode_status.md`** - 🎯 **Recommended starting point**
   - Color-coded overview with 🟢🟠🔴 status indicators
   - Executive summary with statistics
   - Categorized variables by namelist
   - Specific recommendations for improvements

2. **`test_sourcecode_full.md`** - Comprehensive detailed report
   - Complete tables for each namelist
   - All variables with their status across 4 categories
   - Full variable listings and cross-references

## Quick Start

For a quick overview of variable status:
```bash
python3 test_sourcecode.py --report status
```

Then open `test_sourcecode_status.md` to see the color-coded analysis.