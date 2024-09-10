import * as React from 'react';
import { useState, useEffect } from 'react';
import { createRoot } from 'react-dom/client';
import Home from './components/Home.jsx';
import Heatmap from './components/Heatmap.jsx';
import Ordination from './components/Ordination.jsx';

function App() {
    const [currentPage, setCurrentPage] = useState('home');

    useEffect(() => {
        // Listen for the 'navigate-to' IPC message to switch pages
        window.electron.receive('navigate-to', (page) => {
            setCurrentPage(page);
        });
    }, []);

    return (
        <>
            {currentPage == 'home' && <Home />}
            {currentPage == 'heatmap' && <Heatmap />}
            {currentPage == 'ordination' && <Ordination />}
        </>
    );
}

const root = createRoot(document.getElementById('root'));
root.render(<App />);
