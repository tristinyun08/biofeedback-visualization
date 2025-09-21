// FINAL ROBUST VERSION: Works WITH Chrome's security policies.
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
      // A flag to remember if the player was playing before a user-initiated seek.
      let wasPlayingBeforeSeek = false;

      // Event listeners for normal playback.
      audioPlayer.addEventListener('timeupdate', () => updateMarker());
      audioPlayer.addEventListener('play', () => updateMarker());
      audioPlayer.addEventListener('pause', () => updateMarker());

      // After a seek is complete, we check if we need to resume playback.
      audioPlayer.addEventListener('seeked', () => {
        if (wasPlayingBeforeSeek) {
          audioPlayer.play();
        }
        // Always reset the flag and sync the marker after a seek.
        wasPlayingBeforeSeek = false; 
        updateMarker();
      });
      
      updateMarker(); // Initial draw.
      clearInterval(intervalId);

      // --- PLOT -> AUDIO COMMUNICATION (THE SEEKER) ---
      el.on('plotly_click', function(data) {
        if (data && data.points && data.points.length > 0) {
          const clickedTime = data.points[0].x;
          
          if (audioPlayer && typeof clickedTime === 'number' && isFinite(clickedTime)) {
            // 1. Remember if the audio was playing BEFORE we do anything.
            wasPlayingBeforeSeek = !audioPlayer.paused;
            
            // 2. If it was playing, pause it temporarily to prevent glitches.
            if (wasPlayingBeforeSeek) {
              audioPlayer.pause();
            }

            // 3. Set the time. This is a "seek" operation.
            //    This will trigger the 'seeked' event listener when complete.
            audioPlayer.currentTime = clickedTime;
            
            // 4. Manually update the marker for instant visual feedback.
            updateMarker(clickedTime);
          }
        }
      });
    }
  }, 250);
}