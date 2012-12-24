#include <gtk/gtk.h>
#include <glib.h>
#include <stdlib.h>

/*
  gcc gtk_sample2.c `pkg-config --cflags --libs glib-2.0 gtk+-2.0`
*/

char *imagefile = "image2.jpg";

typedef struct rectangle {
  double x;
  double y;
  double width;
  double height;
} Rect;

Rect g_rect = {
  0.0,
  0.0,
  80.0,
  50.0,
};

//Rect *g_accept_points[MAXPOINT];	/* <- using glib.h */
GList* g_reject_points = NULL;

Rect* new_current_g_rect()
{
  Rect *prect = (Rect*)malloc(sizeof(Rect));
  *prect = g_rect;
  return prect;
}

static void
cb_button_accept(GtkWidget *button, gpointer user_data)
{
  printf("%s-%d\n",__FUNCTION__, __LINE__);
  gtk_widget_queue_draw(GTK_WIDGET(user_data));  
  gtk_main_quit();		/* <- exit gtk_main loop */
  g_reject_points = g_list_append(g_reject_points,
				  (gpointer)new_current_g_rect());
}

static void
cb_button_reject(GtkWidget *button, gpointer user_data)
{
  printf("%s-%d\n",__FUNCTION__, __LINE__);
  gtk_widget_queue_draw(GTK_WIDGET(user_data));  
  gtk_main_quit();
}


gboolean cb_expose_event(GtkWidget *widget,
			 GdkEventExpose *event,
			 gpointer user_data)
{
  GdkWindow *drawable = widget->window;
  cairo_t *cr;
  cr = gdk_cairo_create (drawable);

  gdk_cairo_set_source_pixbuf (cr, (GdkPixbuf *)user_data, 0.0, 0.0);
  cairo_paint (cr);
  

  cairo_set_source_rgb (cr, 1.0, 0.0, 0.0);
  /* cairo_rectangle(cairo_t *cr,
     double x,double y, double width,double height); */
  cairo_rectangle (cr, g_rect.x, g_rect.y, g_rect.width, g_rect.height);
  cairo_stroke (cr);
  cairo_destroy (cr);
  return FALSE;
}

int main (int argc, char *argv[])
{
  GtkWidget *window;
  GtkWidget *canvas;
  GdkPixbuf *pixbuf;

  int stride;
  int offset;
  int image_width; 
  int image_height;
  
  gtk_init (&argc, &argv);
  pixbuf = gdk_pixbuf_new_from_file (imagefile,NULL);
  image_width  = gdk_pixbuf_get_width(pixbuf);
  image_height = gdk_pixbuf_get_height(pixbuf);

  window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title (GTK_WINDOW (window),
			"cairo_source_image Sample");
  gtk_widget_set_size_request (window,image_width,image_height*1.2);

  g_signal_connect (G_OBJECT (window), "destroy",
		    G_CALLBACK (gtk_main_quit), NULL);

  canvas = gtk_drawing_area_new ();
  g_signal_connect (G_OBJECT (canvas), "expose-event",
		    G_CALLBACK (cb_expose_event), pixbuf);
  {
    GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
    GtkWidget *button;

    gtk_container_add (GTK_CONTAINER (window), vbox);
    gtk_box_pack_start(GTK_BOX(vbox), canvas, TRUE, TRUE, 0);

    button = gtk_button_new_with_label("accept");
    gtk_widget_set_size_request(GTK_WIDGET(button), 100, 30);
    g_signal_connect (G_OBJECT(button), "clicked",
		      G_CALLBACK(cb_button_accept), canvas);
    gtk_box_pack_start(GTK_BOX(vbox), button,
		       FALSE         , FALSE, 0);

    button = gtk_button_new_with_label("reject");
    gtk_widget_set_size_request(GTK_WIDGET(button), 100, 30);
    g_signal_connect (G_OBJECT(button), "clicked",
		      G_CALLBACK(cb_button_reject), canvas);
    gtk_box_pack_start(GTK_BOX(vbox), button,
		       FALSE        , FALSE, 0);
  }

  gtk_widget_show_all (window);

  for (; g_rect.y < image_height; g_rect.y += g_rect.height){
    for (g_rect.x = 0.0; g_rect.x < image_width; g_rect.x += g_rect.width){
      gtk_main ();
    }
  }

  
  {				/* print get points */
    GList* node = g_reject_points;
    while (node != NULL){
      Rect *data = (Rect *)node->data;
      printf("(%lf,%lf)\n",data->x, data->y);
      node = g_list_next(node);
    }
  }

  g_object_unref (pixbuf);

  return 0;
}
