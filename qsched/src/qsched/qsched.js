import React, { Component } from 'react';
import {Card, CardBody, CardText, Col,
        Form, FormGroup, Label, Input, Button,
        Modal, ModalHeader, ModalBody, ModalFooter,
        Dropdown, DropdownMenu, DropdownItem, DropdownToggle,
        Nav, NavItem, NavLink, Badge} from 'reactstrap'
import classnames from 'classnames'
import queryString from 'query-string'
// to 'download' js objects (csv results)

//import createHistory from 'history/createBrowserHistory'


import { createBrowserHistory } from 'history'

//.createBrowserHistory
import Qlogo from '../Qlogo.png'

let fileDownload = require('js-file-download')

let app = 'https://linuxremote3.bucknell.edu/qsched/qsched'
//let app = 'http://localhost:5000/qsched'
//let app = 'https://3w96mpfgyb.execute-api.us-east-1.amazonaws.com/dev/qsched'

class Qsched extends Component {
  constructor(props){
    super(props)
    this.state = {
      activeTab: 1,
      helpOpen: false,
      helpTopic: "",
      dropdownOpen: false,
      showModal: false,
      errorMessage: "",
      type: 'quant-sin',
      type2: 'quant-poly',
      dims: '',
      jitter2d: '0.7',
      bins: '',
      bias: '1.0',
      evolution: '',
      output_type: '0-start',
      inclusion: true,
      backfill: true,
      linewidth: '1.0',
      linear: '0.7',
      appendcorner: true
    }
    this.help = {
      type: "Is there a ‘one-size-fits-all’ schedule? There are many ways to select and optimize schedules, but if you want an “everyday” schedule that can be used in a variety of settings, experiences suggests that sampling 25-33% of total points, a quant-sin or quant-poly generator, Bias = 1.5 and Evolution = 2.0 generally results in widely applicable schedules.",
      jitter2d: "Jitter (between 0.01-0.99) defines the length  of each edge of a box that is centered in each quantile; thus specifying jitter = 0.7 means that the sample will be jittered in a (-.7 x 0.7) ---> 50% centered box"
    }

    this.samples = {
          1: {
                  'Basic HSQC': {
                    type: 'quant-sin',
                    dims: '512',
                    bins: '128',
                    bias: '1.5',
                    evolution: '3.0',
                    backfill: '15',
                    output_type: '0-start'
                    },
                  'High Res HSQC': {
                    type: 'quant-sin',
                    dims: '1024',
                    bins: '256',
                    bias: '1.5',
                    backfill: '30',
                    evolution: '3.0',
                    output_type: '0-start'}
                },
          2: {
            'Basic HNCA': {
              type: 'quant-sin',
              type2: 'quant-sin',
              dims: '64 32',
              jitter2d: '0.7',
              bins: '32 20',
              bias: '1.5 1.5',
              evolution: '2.0 2.0',
              output_type: '0-start',
              inclusion: true,
              backfill: true,
              appendcorner: true},
            'High Res HNCA': {
              type: 'quant-sin',
              type2: 'quant-sin',
              dims: '128 40',
              jitter2d: '0.7',
              bins: '48 20',
              bias: '1.5 1.5',
              evolution: '2.0 2.0',
              output_type: '0-start',
              inclusion: false,
              backfill: true,
              appendcorner: true
          }
          },
          3: {
            '32 x 128 (25%)': {
              type: 'quant-sin',
              dims: '128',
              bins: '32',
              bias: '1.5',
              evolution:'3.0',
              backfill: '8',
              output_type: '0-start'
            },
            '42 x 128 (33%)':{
              type: 'quant-exp',
              dims: '128',
              bins: '42',
              bias: '1.0',
              evolution: '3.0',
              backfill: '10',
              output_type: '0-start'
            },
            '64 x 256 (25%)':{
              type: 'quant-sin',
              dims: '256',
              bins: '64',
              bias: '1.5',
              evolution: '3.0',
              backfill: '12',
              output_type: '0-start'
            },
            '85 x 256 (33%)':{
              type: 'quant-exp',
              dims: '256',
              bins: '85',
              bias: '1.0',
              evolution: '3.0',
              backfill: '18',
              output_type: '0-start'
            },
            '130 x 512 (25%)':{
              type: 'quant-sin',
              dims: '512',
              bins: '64',
              bias: '1.5',
              evolution: '3.0',
              backfill: '15',
              output_type: '0-start'
            },
            '169 x 512 (33%)':{
              type: 'quant-exp',
              dims: '512',
              bins: '169',
              bias: '1.0',
              evolution: '3.0',
              backfill: '15',
              output_type: '0-start'
            },
            '257 x 1024 (25%)': {
              type: 'quant-sin',
              dims: '1024',
              bins: '256',
              bias: '1.5',
              evolution: '3.0',
              backfill: '30',
              output_type: '0-start'},
            '337 x 1024 (33%)': {
              type: 'quant-exp',
              dims: '1024',
              bins: '337',
              bias: '1.0',
              evolution: '3.0',
              backfill: '30',
              output_type: '0-start'
            }
          }
        }

    this.handleChange = this.handleChange.bind(this)
    this.handleRun = this.handleRun.bind(this)
    this.toggleModal = this.toggleModal.bind(this)
    this.toggleTab = this.toggleTab.bind(this)
    this.checkLocation = this.checkLocation.bind(this)
    //this.history = createHistory()
    this.history = createBrowserHistory()

    this.unlisten = this.history.listen(this.checkLocation)

  }
  checkLocation(loc,act){
      let q = queryString.parse(loc.search)
      q.mode = parseInt(q.mode, 10)
      // console.log('q --- ', q)

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
        let defaultEx = Object.keys(this.samples[1])[0]

        // console.log("bad loc, redirect")
        this.history.replace(
          {
            search:'?mode=1&example=' + encodeURI(defaultEx),
            state: {mode: 1, example: defaultEx}
          })
      }
  }
  componentDidMount(){
    if (this.history.location.search === ""){
      let defaultEx = Object.keys(this.samples[1])[0]
      this.history.replace(
        {
          search:'?mode=1&example=' + encodeURI(defaultEx),
          state: {mode: 1, example: defaultEx}
        })
    }else{
      this.checkLocation(this.history.location)
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

      if (this.state.activeTab == 2){
        // convert boolean to integer
        this.setState({backfill: this.samples[1]['Basic HSQC'].backfill})
      }
      if (tab == 2){
        this.setState({backfill: this.samples[2]['Basic HNCA'].backfill})
      }
      
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
    console.log("change: ", e.target.id, e.target.type)
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
    if (this.state.activeTab === 2){
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
        console.log("got back s", res)

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

    console.log("backfill is", this.state.backfill)
    // console.log('SAMPLES:', samples)
    //<img src="http://www.facstaff.bucknell.edu/drovnyak/QSched_MarkDSR_Logo128W.png" height="100" width="128" alt="Logo"  />
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
            <ModalBody>
              {this.help[this.state.helpTopic]}
            </ModalBody>
            <ModalFooter>
              <Button color="primary" onClick={()=>this.showHelp('')}>OK</Button>
            </ModalFooter>
          </Modal>
          
          <img src={Qlogo} height="100" width="128" alt="Logo"  />

          <Card className="m-4">
            <CardBody>
              <CardText>Welcome to Qsched!
              QSched is extremely flexible and powerful, but it is possible to create schedules that are not appropriate for
              your application.</CardText>
              <CardText> There is no substitute for visually inspecting the schedules you generate.</CardText>
              <CardText></CardText>
              <CardText> For a guide on how to use this site and more, please click <a href="https://sites.google.com/prod/bucknell.edu/drovnyak/softwaredownloads">here.</a></CardText>


              <hr/>
              <Nav tabs>
                <NavItem>
                  <NavLink
                    className={classnames({ active: this.state.activeTab === 1 })}
                    onClick={() => { this.toggleTab(1); }}>
                  1 Dimension
                  </NavLink>
                </NavItem>
                <NavItem>
                  <NavLink
                    className={classnames({ active: this.state.activeTab === 2 })}
                    onClick={() => { this.toggleTab(2); }}>
                  2 Dimensions
                  </NavLink>
                </NavItem>
                <NavItem>
                  <NavLink
                    className={classnames({ active: this.state.activeTab === 3 })}
                    onClick={() => { this.toggleTab(3); }}>
                  One-Click 1D Scheduling
                  </NavLink>
                </NavItem>
              </Nav>
              <hr/>
              {this.state.activeTab !== 3 ?            
              <Dropdown
                id='pop'
                isOpen={this.state.dropdownOpen}
                toggle={()=>this.setState({dropdownOpen:!this.state.dropdownOpen})}>
                <DropdownToggle caret>
                  Example Schedules
                </DropdownToggle>
                <DropdownMenu>
                {samples}
                </DropdownMenu> 
              </Dropdown> :"" }
              {this.state.activeTab === 3 ?
              <Dropdown
                id='pop'
                isOpen={this.state.dropdownOpen}
                toggle={()=>this.setState({dropdownOpen:!this.state.dropdownOpen})}>
                  <DropdownToggle caret>
                    Select One-Click Schedule
                  </DropdownToggle>
                  <DropdownMenu>
                  {samples}  
                  </DropdownMenu>
                </Dropdown> :""}
              <hr/>
              <Form onSubmit={this.handleRun}>
              {this.state.activeTab === 3 ?
              <FormGroup row>
                <Label sm={5}>Type of schedule generator in the 1st dimension</Label>
                <Col sm={6} className="m-auto">
                <Input type="text" id="type" disabled={true}
                  value={this.state.type} onChange={this.handleChange}>
                </Input>
                </Col>
                <Col sm={1} className="m-auto">
                  {'dims' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
              </FormGroup>:""}
              
              {this.state.activeTab === 3 ?
              <FormGroup row>
                <Label sm={5}>Total points comprising Nyquist Grid</Label>
                <Col sm={6} className="m-auto">
                <Input type="text" id="dims" disabled={true}
                  value={this.state.dims} onChange={this.handleChange}>
                </Input>
                </Col>
                <Col sm={1} className="m-auto">
                  {'dims' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
              </FormGroup>:""}
              {this.state.activeTab === 3 ?
              <FormGroup row>
                <Label sm={5}>Number of points to be sampled</Label>
                <Col sm={6} className="m-auto">
                <Input type="text" id="bins" disabled={true}
                  value={this.state.bins} onChange={this.handleChange}>
                </Input>
                </Col>
                <Col sm={1} className="m-auto">
                  {'dims' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
              </FormGroup>:""}

              {this.state.activeTab === 3 ?
              <FormGroup row>
                <Label sm={5}>Bias - affects characteristics of desired PDF</Label>
                <Col sm={6} className="m-auto">
                <Input type="text" id="bias" disabled={true}
                  value={this.state.bias} onChange={this.handleChange}>
                </Input>
                </Col>
                <Col sm={1} className="m-auto">
                  {'dims' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
              </FormGroup>:""}
              {this.state.activeTab === 3 ?
              <FormGroup row>
                <Label sm={5}>Evolution - the degree of T2 that your PDF spans</Label>
                <Col sm={6} className="m-auto">
                <Input type="text" id="evolution" disabled={true}
                  value={this.state.evolution} onChange={this.handleChange}>
                </Input>
                </Col>
                <Col sm={1} className="m-auto">
                  {'dims' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
              </FormGroup>:""}

              {this.state.activeTab !== 3 ?
              <FormGroup row>
                  <Label sm={5}>Type of schedule generator in 1st dimension{" "}</Label>

                  <Col sm={6} className="m-auto">
                    <Input type="select" id="type"
                      value={this.state.type} onChange={this.handleChange}>
                      <option>quant-exp</option>
                      <option>quant-poly</option>
                      <option>quant-sin</option>
                      <option>guass</option>
                      <option>linear</option>
                      <option>noweight</option>
                    </Input>
                    </Col>
                  <Col sm={1} className="m-auto">
                  {'type' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
                </FormGroup> : "" }

                {this.state.activeTab === 2 ?
                <FormGroup row>
                  <Label sm={5}>Type of schedule generator in 2nd dimension</Label>
                  <Col sm={6} className="m-auto">
                    <Input type="select" id="type2"
                      value={this.state.type2} onChange={this.handleChange}>
                      <option>quant-exp</option>
                      <option>quant-poly</option>
                      <option>quant-sin</option>
                      <option>guass</option>
                      <option>linear</option>
                      <option>noweight</option>
                    </Input>
                  </Col>
                  <Col sm={1} className="m-auto">
                  {'type' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
                </FormGroup> : ""
                }
                {this.state.activeTab !== 3 ?
                <FormGroup row>
                  <Label sm={5}>Total points comprising Nyquist grid{" "}</Label>

                      <Col sm={6} className="m-auto">
                    <Input type="text" id="dims"
                      value={this.state.dims} onChange={this.handleChange}/>
                  </Col>
                  <Col sm={1} className="m-auto">
                  {'dims' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
                </FormGroup> : "" }


                  {this.state.activeTab === 2?
                  <FormGroup row>
                  <Label sm={5}>Percent jitter for 2D quantiles</Label>
                  <Col sm={6} className="m-auto">
                    <Input type="text" id="jitter2d"
                      value={this.state.jitter2d} onChange={this.handleChange}/>
                  </Col>
                  <Col sm={1} className="m-auto">
                    {'jitter2d' in this.help ?
                      <Badge color="secondary" onClick={()=>this.showHelp('jitter2d')}> ?</Badge> : ""}
                  </Col>
                </FormGroup> : ""
                }
              {this.state.activeTab !== 3 ?
                <FormGroup row>
                  <Label sm={5}>Number of points to be sampled</Label>
                    <Col sm={6} className="m-auto">
		    <Input type="text" id="bins"
                      value={this.state.bins} onChange={this.handleChange}/>
                  </Col>
                  <Col sm={1} className="m-auto">
                  {'bins' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
                </FormGroup> : "" }

                {this.state.activeTab !== 3 ?
                <FormGroup row>
                  <Label sm={5}>Bias - Affects characteristics of particular PDF</Label>
                   <Col sm={6} className="m-auto">
		    <Input type="text" id="bias"
                      value={this.state.bias} onChange={this.handleChange}/>
                  </Col>
                  <Col sm={1} className="m-auto">
                  {'bias' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
                </FormGroup> : "" }

                {this.state.activeTab !== 3 ?
                <FormGroup row>
                  <Label sm={5}>Evolution - the degree of T2 that your PDF spans</Label>
                    <Col sm={6} className="m-auto">
		    <Input type="text" id="evolution"
                      value={this.state.evolution} onChange={this.handleChange}/>
                  </Col>
                  <Col sm={1} className="m-auto">
                  {'evolution' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
                </FormGroup> : "" }
              {this.state.activeTab === 1 ?
              <FormGroup row>
                <Label sm={5}>Number of points in linear backfill</Label>
                <Col sm={6} className="m-auto">
                  <Input type="text" id="backfill"
                    value={this.state.backfill} onChange={this.handleChange}>
                    </Input>
                </Col>
                <Col sm={1} className="m-auto">
                  {'dims' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
              </FormGroup>:""}
              {this.state.activeTab === 3 ?
              <FormGroup row>
                <Label sm={5}>Number of points in linear backfill</Label>
                <Col sm={6} className="m-auto">
                  <Input type="text" id="backfill" disabled={true}
                    value={this.state.backfill} onChange={this.handleChange}> 
                    </Input>
                </Col>
                <Col sm={1} className="m-auto">
                  {'dims' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
              </FormGroup>:""}
                <FormGroup row>
                  <Label sm={5}>Output type</Label>
		                <Col sm={6} className="m-auto">
                    <Input type="select" id="output_type"
                      value={this.state.output_type} onChange={this.handleChange}>
                      <option>0-start</option>
                      <option>1-start</option>
                    </Input>
                  </Col>
                  <Col sm={1} className="m-auto">
                  {'output_type' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
                </FormGroup>

                {this.state.type === 'guass' ?
                <FormGroup row>
                  <Label sm={5}>Linewidth</Label>
                  <Col sm={6} className="m-auto">
                    <Input type="text" id="linewidth"
                      value={this.state.linewidth} onChange={this.handleChange}/>
                    </Col>
                    <Col sm={1} className="m-auto">
                    {'linewidth' in this.help ?
                      <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                    </Col>
                </FormGroup> : ""}

                {this.state.type === 'linear' ?
                <FormGroup row>
                  <Label sm={5}>Percent of linear sampling</Label>
                  <Col sm={6} className="m-auto">
                    <Input type="text" id="linear"
                      value={this.state.linear} onChange={this.handleChange}/>
                  </Col>
                  <Col sm={1} className="m-auto">
                  {'linear' in this.help ?
                    <Badge color="secondary" onClick={()=>this.showHelp('type')}> ?</Badge> : ""}
                  </Col>
                </FormGroup> : ""}


                {this.state.activeTab===2 ?
                <FormGroup check>
                    <Label check>
                      <Input type="checkbox"  id="inclusion"
                        checked={this.state.inclusion} onChange={this.handleChange}/>
                      Inclusion (edge forcing)
                    </Label>
                  </FormGroup> : ""
                }

                {this.state.activeTab===2 ?
                  <FormGroup check>
                    <Label check>
                      <Input type="checkbox" id="backfill"
                        checked={this.state.backfill} onChange={this.handleChange}/>
                      Backfill
                    </Label>                  </FormGroup> : ""
                }

                {this.state.activeTab===2 ?
                  <FormGroup check>
                    <Label check>
                      <Input type="checkbox" id="appendcorner"
                        checked={this.state.appendcorner} onChange={this.handleChange}/>
                      Append top right corner
                    </Label>
                  </FormGroup> : ""
                }
             
        
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
