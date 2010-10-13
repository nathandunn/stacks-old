var http_req;

function toggle_annotation(id) {
    var div_obj = document.getElementById(id + "_div"); 
    var a_obj   = document.getElementById(id + "_ann"); 
    var frm_obj = document.getElementById(id + "_frm"); 

    // Turn ON the anchor, turn OFF the form
    if (a_obj.style.display == "none") {
        a_obj.style.display = "";
        div_obj.style.display = "none";
        if (a_obj.innerHTML.length > 0) 
            a_obj.innerHTML = "annotate";
    }
    // Turn OFF the anchor, turn ON the form
    else {
        a_obj.style.display   = "none";
        div_obj.style.display = "";
        if (a_obj.innerHTML != "annotate") 
            frm_obj.ext_id.value  = a_obj.innerHTML;
    }
}

function annotate_marker(id) {
    //
    // Fetch the marker annotation
    //
    var form_obj = document.getElementById(id + "_frm"); 
    var url = 
        form_obj.url.value + "&" + 
        "ext_id=" + form_obj.ext_id.value;

    // 
    // Prepare and send XMLHttpRequest Object.
    //
    http_req = false;

    try {
	http_req = new XMLHttpRequest();
    } catch(e) {
	http_req = false;
    }

    if (http_req) {
	http_req.onreadystatechange = process_annotation;
	http_req.open("GET", url, true);
	http_req.send("");
    }

    toggle_annotation(id);
}

function process_annotation() {
    //
    // Possible readyState values:
    // 0 = uninitialized
    // 1 = loading
    // 2 = loaded
    // 3 = interactive
    // 4 = complete
    //
    if (http_req.readyState == 4) {

        // Check that the status is "OK"
        if (http_req.status == 200) {

            var xml_doc = http_req.responseXML;
            var tag_obj = xml_doc.getElementsByTagName("marker_id");
            var txt_obj = xml_doc.getElementsByTagName("text");
            var tag_id  = tag_obj[0].childNodes[0].nodeValue;
            var txt;
            if (txt_obj[0].childNodes.length > 0)
                txt = txt_obj[0].childNodes[0].nodeValue;
            else
                txt = "annotate";

            var a_obj = document.getElementById(tag_id + "_ann");
            a_obj.innerHTML = txt;

        } else {
            alert("There was a problem retrieving the XML data:\n" + http_req.statusText);
        }
    }
}