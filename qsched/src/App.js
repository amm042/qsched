import React, { Component } from 'react';
import {Container} from 'reactstrap';
import Qsched from './qsched/qsched.js';
import Qfooter from './qsched/qfooter.js'
import './App.css';


class App extends Component {
  render() {
    return (
      <div>

      <Container>rt


        <Qsched/>

        <Qfooter/>
      </Container>
      </div>
    );
  }
}

export default App;
