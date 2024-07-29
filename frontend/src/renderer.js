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

// Run plot_relabund_heatmap script
document.getElementById('r-script-form').addEventListener('submit', async (event) => {
    event.preventDefault();
    
    const form = event.target;
    const formData = new FormData(form);
    const params = {};
    
    formData.forEach((value, key) => {
      params[key] = value === 'on' ? true : value; // Convert checkboxes to boolean
    });

   // Add file path to params if available
  if (window.currentFilePath) {
    params.filePath = window.currentFilePath;
  } else {
    console.error('File path is not available');
    return;
  }

  // Add selected ExpObj to params
  const ExpObj = getSelectedExpObj();
  if (ExpObj) {
    params.ExpObj = ExpObj;
  } else {
    console.error('ExpObj is not selected');
  }

  try {
    const result = await window.electron.runRScript(params);
    document.getElementById('output').textContent = result;
  } catch (error) {
    document.getElementById('output').textContent = error;
  }
});
  