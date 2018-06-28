import React, { Component } from 'react';
import {Card, CardBody, CardText, CardTitle,
        Form, FormGroup, Label, Input, Button,
        Modal, ModalHeader, ModalBody, ModalFooter} from 'reactstrap';

let fileDownload = require('js-file-download')

let app = 'http://localhost:5000/qsched'


class Qsched extends Component {
  constructor(props){
    super(props)
    this.state = {
      showModal: false,
      errorMessage: "",
      type: 'quant-sin',
      dims: '1024',
      jitter2d: '0.7',
      bins: '128',
      bias: '1.0',
      evolution: '2.5',
      output_type: 'varian',
      inclusion: true,
      backfill: true,
      linewidth: '1.0',
      linear: '0.7',
      appendcorner: true
    }
    this.handleChange = this.handleChange.bind(this)
    this.handleRun = this.handleRun.bind(this)
    this.toggleModal = this.toggleModal.bind(this)
  }
  toggleModal(){
    this.setState({showModal: !this.state.showModal})
  }
  handleChange(e){
    let k = {}

    k[e.target.id] = e.target.type==='checkbox' ? e.target.checked : e.target.value

    //console.log(k)
    this.setState(k)
  }
  mkcsv(data){
    return data.map(x => {
      if (x.length === undefined)
        return x
      else
        return x.join(", ")
    }).join("\n")
  }
  handleRun(evt){
    console.log(this.state)
    evt.preventDefault()

    fetch(app, {
      body: JSON.stringify(this.state),
      headers: {
        "Content-Type": "application/json"
      },
      method: "POST"
    })
      .then(res => res.json())
      .then(res => {
        console.log("got back s", res)
        if (res.error !== undefined){
          this.setState({
            showModal: true,
            errorMessage: res.error
          })
        }else{
          if (res.length === 2){
            console.log("two files to download")
            fileDownload(this.mkcsv(res[0]), 'schedule.csv')
            fileDownload(this.mkcsv(res[1]), 'schedule_psf.csv')
          }else{
            console.log("one file to download")
            fileDownload(this.mkcsv(res), 'schedule.csv')
          }
        }
      })

  }
  render(){
    return(
        <div>
          <Modal isOpen={this.state.showModal} toggle={this.toggleModal}>
            <ModalHeader toggle={this.toggleModal}>Oops, something went wrong...</ModalHeader>
             <ModalBody>
               {this.state.errorMessage}
             </ModalBody>
             <ModalFooter>
               <Button color="primary" onClick={this.toggleModal}>OK</Button>
             </ModalFooter>
          </Modal>
          <Card className="m-4">
            <CardBody>
              <CardTitle>qsched</CardTitle>
              <CardText>Some intro text goes here. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Donec nec fermentum risus, nec sollicitudin ligula. Aliquam a posuere metus. Sed eget tincidunt nunc. Quisque gravida pellentesque mattis. Cras at quam ullamcorper, placerat magna at, ultricies ex. Proin volutpat orci quis ornare porta. Phasellus non gravida elit. Nunc sed purus et mauris posuere tempus.</CardText>
              <hr/>
              <Form onSubmit={this.handleRun}>
                <FormGroup>
                  <Label>
                    Type of schedule generator
                    <Input type="select" id="type"
                      value={this.state.type} onChange={this.handleChange}>
                      <option>quant-exp</option>
                      <option>quant-poly</option>
                      <option>quant-sin</option>
                      <option>guass</option>
                      <option>linear</option>
                      <option>noweight</option>
                    </Input>
                  </Label>
                </FormGroup>

                <FormGroup>
                  <Label>
                    Indel dimensions
                    <Input type="text" id="dims"
                      value={this.state.dims} onChange={this.handleChange}/>
                  </Label>
                </FormGroup>

                <FormGroup>
                  <Label>
                    Percent jitter for 2D quantiles
                    <Input type="text" id="jitter2d"
                      value={this.state.jitter2d} onChange={this.handleChange}/>
                    <small className="form-text text-muted">(0.7 recommended; between .01-.99)</small>
                  </Label>
                </FormGroup>
                <FormGroup>
                  <Label>
                    Number of points to be sampled along each dimension
                    <Input type="text" id="bins"
                      value={this.state.bins} onChange={this.handleChange}/>
                  </Label>
                </FormGroup>
                <FormGroup>
                  <Label>
                    Bias strength on exponent
                    <Input type="text" id="bias"
                      value={this.state.bias} onChange={this.handleChange}/>
                  </Label>
                </FormGroup>
                <FormGroup>
                  <Label>
                    T2 evolution time for sampling in each indirect dimension
                    <Input type="text" id="evolution"
                      value={this.state.evolution} onChange={this.handleChange}/>
                  </Label>
                </FormGroup>

                <FormGroup>
                  <Label>
                    Output type
                    <Input type="select" id="output_type"
                      value={this.state.output_type} onChange={this.handleChange}>
                      <option>bruker</option>
                      <option>jeol</option>
                      <option>varian</option>
                    </Input>
                  </Label>
                </FormGroup>
                <FormGroup>
                  <label>
                    Linewidth
                    <Input type="text" id="linewidth"
                      value={this.state.linewidth} onChange={this.handleChange}/>
                    <small id="linewidthHelp" className="form-text text-muted">required for guassian option</small>
                  </label>
                </FormGroup>
                <FormGroup>
                  <Label>Percent of linear sampling
                    <Input type="text" id="linear"
                      value={this.state.linear} onChange={this.handleChange}
                      />
                  </Label>
                </FormGroup>

                <FormGroup check>
                  <Label check>
                    <Input type="checkbox"  id="inclusion"
                      checked={this.state.inclusion} onChange={this.handleChange}/>
                    Inclusion (edge forcing)
                  </Label>
                </FormGroup>
                <FormGroup check>
                  <Label check>
                    <Input type="checkbox" id="backfill"
                      checked={this.state.backfill} onChange={this.handleChange}/>
                    Backfill
                  </Label>
                </FormGroup>
                <FormGroup check>
                  <Label check>
                    <Input type="checkbox" id="appendcorner"
                      checked={this.state.appendcorner} onChange={this.handleChange}/>
                    Append top right corner
                  </Label>
                </FormGroup>

              </Form>
              <hr/>

              <Input type="submit" className="btn btn-primary" value="Run"
                onClick={this.handleRun}/>

            </CardBody>
          </Card>
        </div>
    )
  }
}

export default Qsched;
