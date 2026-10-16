# DNA Methylation Analysis Pipeline

Pipeline analisis DNA methylation menggunakan tiga script R untuk analisis utama dan satu script Python untuk membuat laporan HTML.

## File dan Fungsinya

- metadata_template_cust.xlsx: Template metadata yang diberikan kepada customer untuk diisi.
- metadata.csv: Metadata sampel yang digunakan oleh pipeline analisis.
- 01_qc_preprocessing.R: Quality control dan preprocessing data methylation.
- 02_differential_methylation.R: Analisis DMP dan DMR antara kelompok Control dan Case.
- 03_epigenetic_age.R: Estimasi epigenetic age menggunakan Horvath Pan-Tissue Clock.
- methy-ageclock-report.py: Membuat laporan hasil analisis dalam format HTML.
- pandu-logo.png dan pandu-slogan.png: Gambar yang diperlukan untuk laporan HTML.

## Penempatan File

Ketiga script R, metadata.csv, script Python, dan kedua gambar harus disimpan dalam folder kerja yang sama.

Struktur folder:

    folder_analisis/
    ├── 01_qc_preprocessing.R
    ├── 02_differential_methylation.R
    ├── 03_epigenetic_age.R
    ├── methy-ageclock-report.py
    ├── metadata.csv
    ├── pandu-logo.png
    └── pandu-slogan.png

Penting: 
metadata_template.xlsx adalah template untuk customer, bukan file yang dibaca langsung oleh script analisis. Data yang sudah diisi perlu disiapkan sebagai metadata.csv dengan nama kolom yang sesuai.
Ketiga script R menggunakan path file yang ditetapkan di dalam kode. Pastikan idat_dir, metadata_file, dan path input RDS mengarah ke lokasi data yang benar
