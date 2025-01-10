import React from 'react';
import ReactDOM from 'react-dom';
import App from './components/App';

ReactDOM.render(
    <React.StrictMode>
        <App handleNavigateTo={(page) => console.log(`Navigate to ${page}`)} />
    </React.StrictMode>,
    document.getElementById('root')
);

