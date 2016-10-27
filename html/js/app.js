$(document).ready(function() {

  Handlebars.registerHelper('cycle', function(value, index, block) {
    var values = value.split(' ');
    return values[index % values.length];
  });

  Handlebars.registerHelper('json', function(context) {
    return JSON.stringify(context);
  });

  $("#table_id").hide();
  $("form").submit(function(event) {
      event.preventDefault();

      // 09/06/16: I tried cleaning up the code below by doing
      //    ProbAnnoData = new FormData()
      // here, up front, and then using ProbAnnoData.set() to
      // set the values. But that gave a Illegal invocation error
      // on the .ajax call. So I went back to the convoluted
      // code below, lifted from Kris's Kaviar code, apparently
      // for good reason.
     
      // Initialize the input data dictionary, to be
      // used in cgi script as FieldStorage.
      var ProbAnnoData = {
	"template": $("select[name=template]").val()
      };

      var ctType = "application/x-www-form-urlencoded; charset=UTF-8"; 
      var process = "true";
      var taxidCheck = $("input[name=taxid]").val();
      var fileCheck = $("input[name=uploaded_file]").val();

      // if taxid is entered, use that
      if (taxidCheck) {
	ProbAnnoData["taxid"] = taxidCheck;
      }

      // else, if a file is uploaded, use that
      else if (fileCheck) {
	ProbAnnoData["uploaded_file"] = fileCheck;
	var fileInput = $("input[name=uploaded_file]");
	var file = fileInput.get(0).files[0];
	ProbAnnoData = new FormData();
	ProbAnnoData.append('uploaded_file',file);
	ProbAnnoData.append('uploaded_filename',fileCheck);
        // need to reset template since we reset ProbAnnoData
	ProbAnnoData.append('template',$("select[name=template]").val());
	ctType = false;
	process = false;
      }	

      $("#loading").show();
      $.ajax({
	type: "POST",
	url: "http://probanno.systemsbiology.net/cgi-pub/ProbAnno",
	data: ProbAnnoData,
	contentType: ctType,
	processData: process,
	cache: false,
	timeout: 1200000,   // 20 minutes
	//error: function (xhr) {
	   //alert(JSON.stringify(xhr));
	//},
	beforeSend: function(jqXHR, settings) {
	  console.log(settings.url+ '?' + settings.data);
	},
	success: function(data, textStatus, jqXHR){
	  $("#loading").hide();
	  if (data.error) { // if there was an error
	        $("#messages").show();
		if (data.error == '1') {
		     $("#messages").html("<p>ProbAnno algorithm failed</p>");
		} else if (data.error == '2') {
		     $("#messages").html("<p>Could not create output file</p>");
		} else if (data.error == '3') {
		     $("#messages").html("<p>Empty results file</p>");
		} else if (data.error == '7') {
		     $("#messages").html("<p>wget failed (probably invalid Taxonomy ID)</p>");
		} else if (data.error == '8') {
		     $("#messages").html("<p>empty fasta file (probably invalid Taxonomy ID)</p>");
		} else if (data.error == '99') {
		     $("#messages").html("<p>Must specify taxid or file upload</p>");
		}
	  } else if (data.sites.length == 0 ) { // if we have no results
	     $("#messages").show();
	     $("#messages").html("<p>empty ProbAnno results (probably invalid Taxonomy ID)</p>");
	  } else {
	     $("#messages").hide();
		  //display the result in the browser in tsv format
		  // var templateFile = location.origin+"/templates/tsv.html";
		  // $.get(templateFile, function(response) {
		    // var template = response;
		    // var templateScript = Handlebars.compile(template);
		    // var html = templateScript(data);
		    // $("#result").html(html);
		    // $("#result").show();
		  // });

		  // display download link for tsv
		  var downloadTemplateFile = location.origin+"/templates/tsv.download.html";
		  var downloadLinkName = "Download tsv file";
		  $.get(downloadTemplateFile, function(response) {
		    var downloadScript = Handlebars.compile(response);
		    var html = downloadScript(data);
		    // create download link
		    var link = '<a download="probanno.tsv" href="data:text/plain;charset-utf8,'+encodeURIComponent(html)+'">'+downloadLinkName+'</a><br />';
		    $("#downloadLink").html(link);
		  });
		  
		  // display download link for json
		  // downloadTemplateFile = location.origin+"/templates/json.download.html";
		  // var downloadLinkNameJSON = "Download JSON";
		  // $.get(downloadTemplateFile, function(response) {
		    // var downloadScript = Handlebars.compile(response);
		    // var html = downloadScript(data);
		    // // create download link
		    // var link = '<a download="probanno.json" href="data:application/json;charset-utf8,'+encodeURIComponent(html)+'">'+downloadLinkNameJSON+'</a>';
		    // $("#downloadLink").append(link);
		  // });
	      } // end if data
	},
	error: function(jqXHR, textStatus, errorThrown){
	  console.log(jqXHR + "  jqXHR in error");
	  console.log(textStatus + "  textStatus in error");
	  console.log(errorThrown + "  errorThrown in error");
	}
      });
  }); //end form submit
});// end document.ready
