#include <gtk/gtk.h>

/* 
 * gtk_widget_queue_draw(GTK_WIDGET(canvas))
 * で、 on_expose_event 発信する
 * ボタンを押した際にシグナルを結び付けた callback 関数に
 * 再描画したい widget/canvas/drawingarea 等を入れて、
 * gtk_widget_queue_draw(GTK_WIDGET(canvas))
 * として、イベントを発信すれば良いのではないか
 */

GtkWidget *canvas;
GtkWidget *image;
GdkPixbuf *pixbuf;
cairo_surface_t *surface;
char *imagefile = "image.jpg";

typedef struct _rectangle
{
  double x0, y0;
  double x1, y1;
} Rectangle;

Rectangle rect = {50.0, 50.0, 80.0, 80.0};

static void
cb_destroy (GtkWidget *widget, gpointer user_data)
{
  gtk_main_quit ();
}

static void
cb_button_clicked (GtkWidget *button, gpointer user_data)
{
  gtk_main_quit ();
}

static void
cb_button_left (GtkWidget *button, gpointer user_data)
{
  Rectangle *rect = (Rectangle *)user_data;
  rect->x0 -= 10;
  rect->x1 -= 10;

  gtk_widget_queue_draw(GTK_WIDGET(canvas));
  /*gtk_widget_queue_draw(); */
  printf("(%f,%f)\n",rect->x0, rect->x1);
}

static void
cb_button_right (GtkWidget *button, gpointer user_data)
{
  Rectangle *rect = (Rectangle *)user_data;
  rect->x0 += 10;
  rect->x1 += 10;
  //  gtk_widget_queue_draw(canvas);
  gtk_widget_queue_draw(GTK_WIDGET(canvas));
  printf("(%f,%f)\n",rect->x0, rect->x1);
/*   gtk_widget_queue_draw(); */
}

gboolean
cb_expose_event (GtkWidget *widget,
		 GdkEventExpose *event,
		 gpointer user_data)
{
  GdkWindow *drawable = widget->window;
  Rectangle *rect = (Rectangle*)user_data;
  cairo_surface_t *surface = (cairo_surface_t*)user_data;
  cairo_t *cr;
  double x0 = rect->x0;
  double x1 = rect->x1;
  double y0 = rect->y0;
  double y1 = rect->y1;
  
  cr = gdk_cairo_create (drawable);
  
  cairo_set_source_surface(cr, surface, 0.0, 0.0);
  cairo_paint(cr);
  cairo_destroy(cr);

  
  //gdk_cairo_set_source_pixbuf(cr, pixbuf, 0, 0);
  //cairo_paint(cr);
  //cairo_destroy(cr);
  /* cr = gdk_cairo_create (drawable); */
/*   cairo_set_line_width (cr, 10.0); */
/*   cairo_set_source_rgb (cr, 1.0, 0.0, 0.0); */
/*   cairo_rectangle (cr, x0, y0, y0, y1); */
/*   cairo_stroke (cr); */
/*   cairo_destroy (cr); */
  printf("%s\n", __FUNCTION__);
  printf("(%f,%f)\n",rect->x0, rect->x1);

  return FALSE;
}

int main(int argc, char *argv[])
{
  GtkWidget *window;


  gtk_init (&argc, &argv);

  cairo_surface_t *surface = cairo_image_surface_create_from_png(imagefile);

  window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title (GTK_WINDOW (window), "cairo_line_cap Sample");
  //gtk_widget_set_size_request (window, 860, 640);
  gtk_widget_set_size_request (window,
			       cairo_image_surface_get_width(surface)*1.2,
			       cairo_image_surface_get_width(surface)*1.2);
  g_signal_connect(G_OBJECT(window), "destroy",
		   G_CALLBACK(gtk_main_quit), NULL);
  {
    GtkWidget *vbox = gtk_vbox_new (FALSE, 2);
    gtk_container_add (GTK_CONTAINER (window), vbox);
    {
      //      GtkWidget *canvas;
      GtkWidget *hbox;

      canvas = gtk_drawing_area_new();
      g_signal_connect (G_OBJECT (canvas), "expose-event",
			G_CALLBACK (cb_expose_event), surface);

      //image = gtk_image_new_from_file(imagefile);
      //pixbuf = gdk_pixbuf_new_from_file(imagefile, NULL);

      hbox = gtk_hbox_new (FALSE, 2);
      gtk_box_pack_start (GTK_BOX (vbox), canvas, TRUE, TRUE, 0);
      //gtk_box_pack_start (GTK_BOX (vbox), image, TRUE, TRUE, 0);
      gtk_box_pack_start (GTK_BOX (vbox), hbox  , TRUE, TRUE, 0);
      {
	GtkWidget *button;
	button = gtk_button_new_with_label ("left");
	gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, TRUE, 0);
	g_signal_connect (G_OBJECT (button), "clicked",
			  G_CALLBACK (cb_button_left), &rect);

	button = gtk_button_new_with_label ("right");
	gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, TRUE, 0);
	g_signal_connect (G_OBJECT (button), "clicked",
			  G_CALLBACK (cb_button_right), &rect);
      }
    }
  }

  gtk_widget_show_all (window);
  gtk_main();

  return 0;
}
