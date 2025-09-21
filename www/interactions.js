function syncPlotWithAudio(el) {
  let audioPlayer = null;

  const updateMarker = (timeToDraw) => {
    const currentTime = (typeof timeToDraw === 'number' && isFinite(timeToDraw)) 
                        ? timeToDraw 
                        : (audioPlayer ? audioPlayer.currentTime : 0);
    
    if (audioPlayer) {
      const newShape = {
        type: 'line', y0: 0, y1: 1, yref: 'paper',
        x0: currentTime, x1: currentTime,
        line: { color: 'yellow', width: 2 }
      };
      Plotly.relayout(el, { shapes: [newShape] });
    }
  };

  const intervalId = setInterval(() => {
    audioPlayer = document.getElementById('audioplayer');
    if (audioPlayer) {
      let wasPlayingBeforeSeek = false;

      audioPlayer.addEventListener('timeupdate', () => updateMarker());
      audioPlayer.addEventListener('play', () => updateMarker());
      audioPlayer.addEventListener('pause', () => updateMarker());

      audioPlayer.addEventListener('seeked', () => {
        if (wasPlayingBeforeSeek) {
          audioPlayer.play();
        }
        wasPlayingBeforeSeek = false; 
        updateMarker();
      });
      
      updateMarker();
      clearInterval(intervalId);

      el.on('plotly_click', function(data) {
        if (data && data.points && data.points.length > 0) {
          const clickedTime = data.points[0].x;
          
          if (audioPlayer && typeof clickedTime === 'number' && isFinite(clickedTime)) {
            wasPlayingBeforeSeek = !audioPlayer.paused;
            
            if (wasPlayingBeforeSeek) {
              audioPlayer.pause();
            }

            audioPlayer.currentTime = clickedTime;
            updateMarker(clickedTime);
          }
        }
      });
    }
  }, 250);
}