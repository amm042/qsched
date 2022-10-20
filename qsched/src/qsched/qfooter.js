import React, { Component } from 'react';
// import {Container,Row, Col} from 'reactstrap'
//import {Row} from 'reactstrap'
class QFooter extends Component{
  render() {
    return(

      <footer className="footer">
        <div className="h-100 m-auto">
          <p>Copyright 2017, 2018 D. Levi Craft; David Rovnyak, Virginia G. Rovnyak.</p>
          <p>Web implementation by <a href="https://www.eg.bucknell.edu/~amm042/">Prof. Alan Marchiori</a>, <a href="https://www.linkedin.com/in/lucille-cullen-083698251/">Lucy E.Cullen</a> and <a href="https://www.linkedin.com/in/mark-roginkin-915259131/">Mark Roginkin</a>.</p>
          <p>If you use QSched please cite: D. Levi Craft, Reilly E. Sonstrom, Virginia G. Rovnyak, David Rovnyak, “Nonuniform sampling by quantiles”, J. Magn. Reson. 288, 109-121 (2018).</p>
        </div>
      </footer>

    )
  }
}

export default QFooter;
