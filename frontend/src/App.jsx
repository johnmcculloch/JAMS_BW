import React from 'react';
import './renderer/shared/styles/index.css';

function App() {
  return (
    <div>
      <h1>Welcome to the JAMSBeta Analysis App</h1>
      <button id="ordination-btn">Go to Ordination Plot Analysis</button>
      <button id="heatmap-btn">Go to Relabund Heatmap Analysis</button>
      <button id="alphadiversity-btn">Go to Alpha Diversity Plot Analysis</button>
      <button id="relabundfeatures-btn">Go to Relabund Features Analysis</button>
    </div>
  );
}

export default App;
