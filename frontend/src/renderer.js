document.getElementById('rdata-file').addEventListener('change', async (event) => {
  const file = event.target.files[0];
  if (file) {
    const filePath = file.path;
    try {
      const objects = await window.electron.loadRDataFile(filePath);
      const dropdown = document.getElementById('summarized-object');
      dropdown.innerHTML = ''; // Clear previous options
      objects.forEach(obj => {
        const option = document.createElement('option');
        option.value = obj;
        option.textContent = obj;
        dropdown.appendChild(option);
      });
    } catch (error) {
      console.error('Failed to load .RData file:', error);
    }
  }
});


document.getElementById('r-script-form').addEventListener('submit', async (event) => {
    event.preventDefault();
    
    const form = event.target;
    const formData = new FormData(form);
    const params = {};
    
    formData.forEach((value, key) => {
      params[key] = value === 'on' ? true : value; // Convert checkboxes to boolean
    });
  
    try {
      const result = await window.electron.runRScript(params);
      document.getElementById('output').textContent = result;
    } catch (error) {
      document.getElementById('output').textContent = error;
    }
  });
  