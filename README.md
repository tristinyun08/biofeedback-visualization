# Biofeedback Visualization

A comprehensive R Shiny application for analyzing the relationship between physiological data (heart rate variability) and audio performance data. This tool provides real-time visualization and analysis of how musical performance affects autonomic nervous system responses.

## Live Demo

**[Try the application online](https://tristin.shinyapps.io/Biofeedback_Visualization/)**

The application is hosted on ShinyApps.io and ready to use with your own HRV and audio data.

## Features

### Audio Analysis
- **Waveform Visualization**: Real-time audio waveform display with interactive seeking
- **Musical Features**: Analysis of loudness (RMS), pitch (fundamental frequency), and brightness (spectral centroid)
- **Audio Synchronization**: Click on plots to seek through audio playback

### Physiological Analysis
- **Heart Rate Variability (HRV)**: Rolling RMSSD calculation with artifact correction
- **Heart Rate Monitoring**: Real-time BPM tracking with min/max/average summaries
- **Stress Zone Analysis**: Relative stress/recovery zones based on session data

### Interactive Visualizations
- **Combined Performance View**: Synchronized audio waveform and heart rate plots
- **Correlation Analysis**: Cross-correlation between audio features and physiological metrics
- **Joint Analysis Plots**: Dual-axis visualizations showing relationships between metrics
- **Dynamic Correlation Table**: Adjustable lag analysis with real-time correlation coefficients

## Installation

### Prerequisites
- R (version 4.0 or higher)
- RStudio (recommended)

### Required R Packages
```r
install.packages(c(
  "shiny",
  "shinyjs", 
  "bslib",
  "ggplot2",
  "dplyr",
  "tuneR",
  "seewave",
  "av",
  "plotly",
  "htmlwidgets",
  "shinycssloaders"
))
```

### Running the Application
1. Clone this repository
2. Open `app.R` in RStudio
3. Click "Run App" or execute:
```r
shiny::runApp()
```

## Usage

### Data Requirements
- **HRV Data**: Text file containing RR intervals in milliseconds (one per line)
- **Audio Data**: Audio file in supported formats (.m4a, .wav, .mp3, etc.)

### Analysis Workflow
1. Upload your HRV data file (.txt format)
2. Upload your performance audio file
3. Set the analysis duration (10-300 seconds)
4. Click "Run Analysis" to process the data
5. Explore the interactive visualizations across multiple tabs

### Visualization Tabs

#### Performance (Audio + HR)
- Synchronized audio waveform and heart rate visualization
- Summary statistics (min, max, average heart rate)
- Interactive seeking through audio playback

#### Relative Stress Sync
- RMSSD visualization with relative stress zones
- Session-based stress/recovery classification
- Dynamic threshold calculation

#### Audio Features Sync
- Loudness (RMS energy) over time
- Sound brightness (spectral centroid)
- Mean pitch (fundamental frequency)

#### Joint Analysis & Correlation
- Cross-correlation analysis with adjustable lag
- Dual-axis plots showing physiological and audio metrics
- Real-time correlation coefficient calculations

## Technical Details

### Audio Processing
- Automatic format conversion using the `av` package
- 1-second window analysis for feature extraction
- RMS energy calculation for loudness
- Fundamental frequency detection for pitch analysis
- Spectral centroid calculation for brightness

### HRV Analysis
- Artifact detection and correction (20% threshold)
- Rolling RMSSD calculation with 30-second windows
- Instantaneous heart rate calculation
- Time-series binning for synchronization

### Interactive Features
- Plotly-based interactive visualizations
- Audio-player synchronization
- Click-to-seek functionality
- Real-time correlation analysis

## File Structure

```
biofeedback-visualization/
├── appv2.R              # Main Shiny application
├── www/
│   ├── styles.css       # Custom styling
│   ├── interactions.js  # Audio-plot synchronization
│   └── current_audio.wav # Temporary audio file
├── rsconnect/           # Deployment configuration
└── README.md           # This file
```

## Acknowledgments

- Built with R Shiny framework
- Audio processing powered by `tuneR` and `seewave` packages
- Interactive visualizations using `plotly`
- Styling with `bslib` and custom CSS


