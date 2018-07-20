import React, { Component } from 'react';
import {Card, CardBody, CardText, CardTitle,
        Form, FormGroup, Label, Input, Button,
        Modal, ModalHeader, ModalBody, ModalFooter,
        Dropdown, DropdownMenu, DropdownItem, DropdownToggle,
        Nav, NavItem, NavLink, Badge } from 'reactstrap'
import classnames from 'classnames'
import queryString from 'query-string'
// to 'download' js objects (csv results)

import createHistory from 'history/createBrowserHistory'

let fileDownload = require('js-file-download')

//let app = 'http://localhost:5000/qsched'
let app = 'https://3w96mpfgyb.execute-api.us-east-1.amazonaws.com/dev/qsched'

class Qsched extends Component {
  constructor(props){
    super(props)
    this.state = {
      activeTab: '1',
      helpOpen: false,
      helpTopic: "",
      dropdownOpen: false,
      showModal: false,
      errorMessage: "",
      type: 'quant-sin',
      type2: 'quant-poly',
      dims: '',
      jitter2d: '',
      bins: '',
      bias: '1.0',
      evolution: '',
      output_type: 'varian',
      inclusion: true,
      backfill: true,
      linewidth: '1.0',
      linear: '0.7',
      appendcorner: true
    }
    this.help = {
      type: "Some help about type.",
      dims: "Some help about dims."
    }
    this.samples = {
      '1': {
              'HSQC': {
                type: 'quant-sin',
                dims: '1024',
                bins: '128',
                bias: '1.0',
                evolution: '2.5',
                output_type: 'varian',
                inclusion: true,
                backfill: true,
                appendcorner: true},
              'HSQC num 2': {
                type: 'quant-sin',
                dims: '55',
                bins: '44',
                bias: '1.0',
                evolution: '2.5',
                output_type: 'varian',
                inclusion: true,
                backfill: true,
                appendcorner: true}
            },
      '2': {
        'TOCSY': {
          type: 'quant-poly',
          type2: 'quant-poly',
          dims: '90 40',
          jitter2d: '0.8',
          bins: '24 12',
          bias: '1.5 1.5',
          evolution: '3.14 1.0',
          output_type: 'varian',
          inclusion: false,
          backfill: true,
          appendcorner: true
        }
      }
    }
    this.handleChange = this.handleChange.bind(this)
    this.handleRun = this.handleRun.bind(this)
    this.toggleModal = this.toggleModal.bind(this)
    this.toggleTab = this.toggleTab.bind(this)
    this.history = createHistory()

    this.unlisten = this.history.listen((loc, act)=> {
      console.log('loc', loc)
      let q = queryString.parse(loc.search)

      if (('example' in q) &&
          ('mode' in q) &&
          (q.mode in this.samples) &&
          (q.example in this.samples[q.mode])){

        let s = Object.assign({activeTab:q.mode},
          this.samples[q.mode][q.example])

        this.setState(s)
      }else{
        // url was invalid, go to the default
        //document.location.search = '?mode=1&example=HSQC'
        console.log("bad loc, redirect")
        this.history.replace(
          {
            search:'?mode=1&example=HSQC',
            state: {mode: 1, example: 'HSQC'}
          })
      }
    })

  }
  componentDidMount(){
    if (this.history.location.search === ""){
      this.history.replace(
        {
          search:'?mode=1&example=HSQC',
          state: {mode: 1, example: 'HSQC'}
        })
    }else{

    }

  }
  showHelp(topic){
    if (topic in this.help){
      this.setState({
        helpOpen: true,
        helpTopic: topic
      })
    }else{
      this.setState({helpOpen:false})
    }
  }
  toggleTab(tab) {
    if (this.state.activeTab !== tab) {
      this.setState({
        activeTab: tab
      });
    }
  }
  toggleModal(){
    this.setState({showModal: !this.state.showModal})
  }
  handleChange(e){
    let k = {}
    k[e.target.id] = e.target.type==='checkbox' ? e.target.checked : e.target.value
    this.setState(k)
  }
  mkcsv(data){
    // returns a CSV-like string from a JS array object
    return data.map(x => {
      if (x.length === undefined)
        return x
      else
        return x.join(", ")
    }).join("\n")
  }
  handleRun(evt){
    //console.log(this.state)
    evt.preventDefault()

    // fixup for 2d schedule generator function
    let q = Object.assign({}, this.state)
    if (this.state.activeTab === '2'){
      q.type = q.type + ' ' + q.type2
    }
    // POST the arguments to the server and wait for a response
    fetch(app, {
      body: JSON.stringify(q),
      headers: {
        "Content-Type": "application/json"
      },
      method: "POST"
    })
      .then(res => res.json())
      .then(res => {
        //console.log("got back s", res)

        if (res.error !== undefined){
          // something went wrong, show the modal dialog with the message
          this.setState({
            showModal: true,
            errorMessage: res.error
          })
        }else{
          // got the result (should set success:true or something in the future)
          // to be sure, but we didn't get an error at least!
          // convert the JS object to a CSV string and trigger a download
          // from the browser.
          if (res.length === 2){
            //console.log("two files to download")
            fileDownload(this.mkcsv(res[0]), 'schedule.csv')
            fileDownload(this.mkcsv(res[1]), 'schedule_psf.csv')
          }else{
            //console.log("one file to download")
            fileDownload(this.mkcsv(res), 'schedule.csv')
          }
        }
      })

  }
  render(){
    // returns a form containing all of the possible input arguments to
    // the scheduler.
    let samples= Object.keys(this.samples[this.state.activeTab]).map((x)=>{
      console.log("smap--", x)


      return <DropdownItem key={x} tag="a"
        onClick={()=>{this.history.push(
          {
            search: '?mode='+this.state.activeTab+'&example='+x,
            state: {
              mode: this.state.activeTab,
              example: x
            }
          })}}>{x}</DropdownItem>
    })

    console.log("SAMPLES:", samples);
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
          <Modal isOpen={this.state.helpOpen}>
            <ModalHeader>
              {this.state.helpTopic}
            </ModalHeader>
            <ModalBody>
              {this.help[this.state.helpTopic]}
            </ModalBody>
            <ModalFooter>
              <Button color="primary" onClick={()=>this.showHelp('')}>OK</Button>
            </ModalFooter>
          </Modal>

          <Card className="m-4">
            <CardBody>
              <CardTitle>qsched</CardTitle>
              <CardText>This is a public beta release of a NUS scheduling routine.</CardText>
              <hr/>
              <Nav tabs>
                <NavItem>
                  <NavLink
                    className={classnames({ active: this.state.activeTab === '1' })}
                    onClick={() => { this.toggleTab('1'); }}>
                  1 Dimension
                  </NavLink>
                </NavItem>
                <NavItem>
                  <NavLink
                    className={classnames({ active: this.state.activeTab === '2' })}
                    onClick={() => { this.toggleTab('2'); }}>
                  2 Dimensions
                  </NavLink>
                </NavItem>
              </Nav>
              {this.state.helpOpen ? this.help[this.state.helpTopic] : ""}
              <hr/>
              <Dropdown
                id='pop'
                isOpen={this.state.dropdownOpen}
                toggle={()=>this.setState({dropdownOpen:!this.state.dropdownOpen})}>
                <DropdownToggle caret>
                  Examples
                </DropdownToggle>
                <DropdownMenu>
                  {samples}
                </DropdownMenu>

              </Dropdown>
              <hr/>
              <Form onSubmit={this.handleRun}>
                <FormGroup>
                  <Label>
                    Type of schedule generator{" "}
                    {'type' in this.help ?
                      <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
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
                {this.state.activeTab==='2' ?
                <FormGroup>
                  <Label>
                    2nd Type of schedule generator
                    <Input type="select" id="type2"
                      value={this.state.type2} onChange={this.handleChange}>
                      <option>quant-exp</option>
                      <option>quant-poly</option>
                      <option>quant-sin</option>
                      <option>guass</option>
                      <option>linear</option>
                      <option>noweight</option>
                    </Input>
                  </Label>
                </FormGroup> : ""
                }
                <FormGroup>
                  <Label>
                    Indel dimensions{" "}
                    {'dims' in this.help ?
                      <Badge color="secondary" onClick={()=>this.showHelp('dims')}> ?</Badge> : ""}
                    <Input type="text" id="dims"
                      value={this.state.dims} onChange={this.handleChange}/>
                  </Label>
                </FormGroup>

                <FormGroup>
                  {this.state.activeTab==='2'?
                  <Label>
                    Percent jitter for 2D quantiles
                    <Input type="text" id="jitter2d"
                      value={this.state.jitter2d} onChange={this.handleChange}/>
                    <small className="form-text text-muted">(0.7 recommended; between .01-.99)</small>
                  </Label> : ""}
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

                {this.state.type === 'linear' ?
                <FormGroup>
                  <label>
                    Linewidth
                    <Input type="text" id="linewidth"
                      value={this.state.linewidth} onChange={this.handleChange}/>
                    <small id="linewidthHelp" className="form-text text-muted">required for guassian option</small>
                  </label>
                </FormGroup> : ""}

                {this.state.type === 'linear' ?
                <FormGroup>
                  <Label>Percent of linear sampling
                    <Input type="text" id="linear"
                      value={this.state.linear} onChange={this.handleChange}
                      />
                  </Label>
                </FormGroup> : ""}

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
