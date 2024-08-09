// ordination.js

document.getElementById('home-btn').addEventListener('click', () => {
    // send IPC message to main process to navigate to heatmap page
    window.electron.send('navigate-to', 'renderer/home/home.html');
});

//document.getElementById('openFileLocation').addEventListener('click', () => {
  //  window.electron.send('open-file-location');
  //});

  // Listen for the param-str event and display paramStr
window.electron.onParamStr((paramStr) => {
    const paramStrOutput = document.getElementById('param-str-output');
    paramStrOutput.textContent = `Generated Command: ${paramStr}`;
  });
  
  // Loading R Data File for summarized experiment object selection
  document.getElementById('rdata-file').addEventListener('change', async (event) => {
    const file = event.target.files[0];
    if (file) {
      const filePath = file.path;
      try {
        // Pass file path directly to function
        const objects = await window.electron.loadRDataFile(filePath);
        console.log("Loaded objects:", objects)
        const dropdown = document.getElementById('ExpObj');
        if (dropdown) {
          dropdown.innerHTML = ''; // Clear previous options
          objects.forEach(obj => {
            const option = document.createElement('option');
            option.value = obj;
            option.textContent = obj;
            dropdown.appendChild(option);
          });
        } else {
          console.error('Dropdown element with id "ExpObj" not found');
        }
        // Store filePath in a global variable if needed later
        window.currentFilePath = filePath;
      } catch (error) {
        console.error('Failed to load .RData file:', error);
      }
    }
  });


    // Function to get the selected ExpObj
    function getSelectedExpObj() {
        const dropdown = document.getElementById('ExpObj');
        return dropdown ? dropdown.value : null;
      }
