#!/usr/bin/env python3
"""
uDALES Variable Analysis - Generate All Reports

This script generates all three types of reports:
1. Full detailed report with all variables per namelist
2. Compact summary report 
3. Color-coded status report with warnings and errors (recommended)
"""

import subprocess
import sys
from pathlib import Path


def main():
    """Generate all reports in sequence."""
    
    print("🔍 uDALES Variable Analysis - Generating All Reports")
    print("=" * 60)
    
    script_dir = Path(__file__).parent
    
    # List of scripts to run
    scripts = [
        ("generate_reports.py", "Full & Compact Reports"),
        ("generate_status_report.py", "Status Report (Color-coded)")
    ]
    
    for script, description in scripts:
        script_path = script_dir / script
        if not script_path.exists():
            print(f"❌ Error: {script} not found")
            continue
            
        print(f"\n📊 Generating {description}...")
        print("-" * 40)
        
        try:
            result = subprocess.run([sys.executable, str(script_path)], 
                                  cwd=script_dir, 
                                  capture_output=True, 
                                  text=True)
            
            if result.returncode == 0:
                print(result.stdout)
                print(f"✅ {description} completed successfully")
            else:
                print(f"❌ Error generating {description}:")
                print(result.stderr)
                
        except Exception as e:
            print(f"❌ Exception running {script}: {e}")
    
    print("\n" + "=" * 60)
    print("📋 Generated Reports Summary:")
    print("=" * 60)
    
    # List generated files
    report_files = [
        ("variable_index_full.md", "Full Report", "Detailed tables per namelist"),
        ("variable_index_compact.md", "Compact Report", "Executive summary"),
        ("variable_status_report.md", "Status Report", "🟢🟠🔴 Color-coded with warnings/errors")
    ]
    
    for filename, title, description in report_files:
        filepath = script_dir / filename
        if filepath.exists():
            size_kb = filepath.stat().st_size / 1024
            print(f"✅ {title:15} | {filename:25} | {size_kb:6.1f} KB | {description}")
        else:
            print(f"❌ {title:15} | {filename:25} | Missing")
    
    print("\n🎯 Recommended: Start with 'variable_status_report.md' for color-coded overview")
    print("📖 Full documentation: VARIABLE_INDEXER_README.md")


if __name__ == "__main__":
    main()