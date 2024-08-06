// heatmap.js


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
      const element = form.elements[key]; // get form element
      if (element.type === 'checkbox') {
        params[key] = element.checked; // TRUE if checked, false if unchecked
      } else if (key === 'colcategories') {
        // Split the input string by commas and trim whitespace
        const colcategoriesArray = value.split(',').map(item => item.trim());
        // Convert the array to R vector format
        params[key] = `c(${colcategoriesArray.map(item => `"${item}"`).join(', ')})`;
      } else if (value === 'TRUE' || value === 'FALSE') {
        params[key] = (value === 'TRUE'); // Convert TRUE/FALSE to boolean
      } else if (!isNaN(value) && value !== '') {
        params[key] = parseInt(value); // Convert number strings to integers
      } else {
        params[key] = value === '' ? null : value; // convert empty strings to null
      }
    });
  
    // Add unchecked boxes as false
    const checkboxes = form.querySelectorAll('input[type="checkbox"]');
    checkboxes.forEach((checkbox) => {
      if (!checkbox.checked && !(checkbox.name in params)){
        params[checkbox.name] = false;
      }
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
      return;
    }
  
    try {
      const result = await window.electron.runRScript(params);
      document.getElementById('output').textContent = result.stdout;
  
  // Display the generated image or PDF
  const heatmapImage = document.getElementById('heatmap-image');
  const fileType = result.imagePath.split('.').pop().toLowerCase(); // Get file extension
  if (fileType === 'pdf') {
    await renderPDF(result.imagePath);
  } else {
    heatmapImage.src = `file://${result.imagePath}`;
    heatmapImage.style.display = 'block';
    heatmapImage.onload = () => {
      document.body.scrollTop = document.body.scrollHeight;
    };
  }
  } catch (error) {
  document.getElementById('output').textContent = error;
  document.getElementById('heatmap-image').style.display = 'none';
  }
  });
  
  // Function to render a PDF
  async function renderPDF(pdfPath) {
  const pdfContainer = document.getElementById('pdf-container');
  pdfContainer.innerHTML = ''; // Clear previous content
  
  try {
  const pdf = await pdfjsLib.getDocument(pdfPath).promise;
  const numPages = pdf.numPages;
  
  for (let pageNum = 1; pageNum <= numPages; pageNum++) {
    const page = await pdf.getPage(pageNum);
    const canvas = document.createElement('canvas');
    const context = canvas.getContext('2d');
    const viewport = page.getViewport({ scale: 1.5 });
    canvas.height = viewport.height;
    canvas.width = viewport.width;
  
    pdfContainer.appendChild(canvas);
  
    await page.render({ canvasContext: context, viewport: viewport }).promise;
  }
  } catch (error) {
  console.error('Error rendering PDF:', error);
  }
  }
  
  
     