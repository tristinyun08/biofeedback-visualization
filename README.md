# Biofeedback Visualization

A comprehensive R Shiny application for analyzing the relationship between physiological data (heart rate variability) and audio performance data. This tool provides real-time visualization and analysis of how musical performance affects autonomic nervous system responses.

## Features

### 🎵 Audio Analysis
- **Waveform Visualization**: Real-time audio waveform display with interactive seeking
- **Musical Features**: Analysis of loudness (RMS), pitch (fundamental frequency), and brightness (spectral centroid)
- **Audio Synchronization**: Click on plots to seek through audio playback

### ❤️ Physiological Analysis
- **Heart Rate Variability (HRV)**: Rolling RMSSD calculation with artifact correction
- **Heart Rate Monitoring**: Real-time BPM tracking with min/max/average summaries
- **Stress Zone Analysis**: Relative stress/recovery zones based on session data

### 📊 Interactive Visualizations
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
2. Open `appv2.R` in RStudio
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

## Deployment

The application can be deployed to:
- **ShinyApps.io**: Use the `rsconnect` folder configuration
- **RStudio Connect**: Direct deployment from RStudio
- **Local server**: Run locally for development and testing

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Built with R Shiny framework
- Audio processing powered by `tuneR` and `seewave` packages
- Interactive visualizations using `plotly`
- Styling with `bslib` and custom CSS

## Support

For issues and questions, please open an issue on GitHub or contact the development team.
